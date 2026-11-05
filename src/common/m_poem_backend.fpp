!> CPU or NVIDIA OpenACC/double bridge to the unchanged POEM model.
#:include 'macros.fpp'
module m_poem_backend

    use iso_c_binding, only: c_int, c_double, c_char, c_null_char
    use iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use m_precision_select, only: wp
#ifdef MFC_OpenACC
    use openacc
#endif

    implicit none

    private
    public :: s_poem_initialize, s_poem_decode, s_poem_advance, s_poem_advance_active, poem_Aref
    public :: s_poem_decode_batch, s_poem_decode_cells, s_poem_write_profile, s_poem_invalidate_cells
    public :: s_poem_decode_device, s_poem_advance_device, poem_dimension
    logical  :: initialized = .false.
    integer  :: poem_dimension = 0
    real(wp) :: poem_Aref
    $:GPU_DECLARE(create='[poem_Aref, poem_dimension]')
    logical                     :: profile_enabled = .false.
    integer(int64)              :: decode_calls = 0, decode_states = 0, source_calls = 0, source_states = 0
    integer(int64)              :: cell_cache_hits = 0, cell_cache_misses = 0, cell_padding_states = 0
    real(wp)                    :: decode_seconds = 0._wp, source_seconds = 0._wp
    integer(int64), allocatable :: cell_bits(:,:)
    real(wp), allocatable       :: cell_moments(:,:)
    logical, allocatable        :: cell_valid(:)

    interface
        integer(c_int) function c_state_dimension() bind(C, name='mfc_poem_state_dimension')
            import c_int

        end function c_state_dimension
        integer(c_int) function c_init_device(root, library, device) bind(C, name='mfc_poem_init_device')
            import c_int, c_char
            character(c_char), intent(in) :: root(*), library(*)
            integer(c_int), value         :: device

        end function c_init_device
        integer(c_int) function c_decode_device(count, y, moments, stream) bind(C, name='mfc_poem_decode_device')
            import c_int, c_double, int64
            integer(c_int), value       :: count
            real(c_double), intent(in)  :: y(*)
            real(c_double), intent(out) :: moments(*)
            integer(int64), value       :: stream

        end function c_decode_device
        integer(c_int) function c_step_device(count, active, pressure, h, y, following, ratio, root, stream) bind(C, &
                & name='mfc_poem_step_device')
            import c_int, c_double, int64
            integer(c_int), value       :: count
            integer(c_int), intent(in)  :: active(*)
            real(c_double), value       :: h
            real(c_double), intent(in)  :: pressure(*), y(*)
            real(c_double), intent(out) :: following(*), ratio(*), root(*)
            integer(int64), value       :: stream

        end function c_step_device
        integer(c_int) function c_init(root, library) bind(C, name='mfc_poem_init')
            import c_int, c_char
            character(c_char), intent(in) :: root(*), library(*)

        end function c_init
        integer(c_int) function c_constants(values) bind(C, name='mfc_poem_constants')
            import c_int, c_double
            real(c_double), intent(out) :: values(*)

        end function c_constants
        integer(c_int) function c_decode(count, y, moments, valid) bind(C, name='mfc_poem_decode')
            import c_int, c_double
            integer(c_int), value       :: count
            real(c_double), intent(in)  :: y(*)
            real(c_double), intent(out) :: moments(*)
            integer(c_int), intent(out) :: valid(*)

        end function c_decode
        integer(c_int) function c_step(count, pressure, h, y, following, moments, valid, ratio, root) bind(C, name='mfc_poem_step')
            import c_int, c_double
            integer(c_int), value       :: count
            real(c_double), value       :: h
            real(c_double), intent(in)  :: pressure(*), y(*)
            real(c_double), intent(out) :: following(*), moments(*), ratio(*), root(*)
            integer(c_int), intent(out) :: valid(*)

        end function c_step
        integer(c_int) function c_error(message, length) bind(C, name='mfc_poem_error')
            import c_int, c_char
            integer(c_int), value          :: length
            character(c_char), intent(out) :: message(*)

        end function c_error
        subroutine s_poem_abort(message)

            character(len=*), intent(in) :: message

        end subroutine s_poem_abort
    end interface

contains
    subroutine s_poem_check(code)

        integer(c_int), intent(in)       :: code
        integer(c_int)                   :: ignored
        integer                          :: i
        character(kind=c_char, len=2048) :: message

        if (code == 0) return
        message = ' '
        ignored = c_error(message, int(len(message), c_int))
        do i = 1, len(message)
            if (message(i:i) == c_null_char) message(i:) = ' '
        end do
        call s_poem_abort('POEM backend: ' // trim(message))

    end subroutine s_poem_check

    subroutine s_poem_initialize()

        character(len=2048) :: root, library
        character(len=16)   :: profile_setting
        real(c_double)      :: values(6)

        if (initialized) return
        if (storage_size(0._wp) /= 64 .or. storage_size(0._c_double) /= 64) &
            & call s_poem_abort('The batched POEM bridge requires double precision.')
        call s_poem_invalidate_cells()
        call get_environment_variable('POEM_MFC_PROFILE', profile_setting)
        profile_enabled = trim(profile_setting) == 'T'
        call get_environment_variable('POEM_REPO_ROOT', root)
        call get_environment_variable('POEM_MFC_LIBRARY', library)
        if (len_trim(root) == 0 .or. len_trim(library) == 0) &
            & call s_poem_abort('Set POEM_REPO_ROOT and POEM_MFC_LIBRARY for the CPU POEM backend.')
#ifdef MFC_OpenACC
        call s_poem_check(c_init_device(trim(root) // c_null_char, trim(library) // c_null_char, &
                          & int(acc_get_device_num(acc_device_nvidia), c_int)))
#else
        call s_poem_check(c_init(trim(root) // c_null_char, trim(library) // c_null_char))
#endif
        poem_dimension = c_state_dimension()
        if (poem_dimension < 4 .or. mod(poem_dimension, &
            & 2) /= 0) call s_poem_abort('POEM backend returned an invalid state dimension.')
        call s_poem_check(c_constants(values))
        poem_Aref = values(1)
        $:GPU_UPDATE(device='[poem_Aref, poem_dimension]')
        initialized = .true.

    end subroutine s_poem_initialize

    !> All arrays must already be present on the device. Synchronous boundary deliberately precedes stream-0 handoff; no bulk host
    !! transfers occur.
    subroutine s_poem_decode_device(count, y, moments)

        integer, intent(in)   :: count
        real(wp), intent(in)  :: y(poem_dimension, count)
        real(wp), intent(out) :: moments(4, count)
        integer(c_int)        :: code
        integer(int64)        :: started

        call s_poem_initialize()
#ifdef MFC_OpenACC
        !$acc wait
        decode_calls = decode_calls + 1_int64
        decode_states = decode_states + int(count, int64)
        if (profile_enabled) call system_clock(started)
        #:call GPU_HOST_DATA(use_device_addr='[y, moments]')
            code = c_decode_device(int(count, c_int), y, moments, 0_int64)
        #:endcall GPU_HOST_DATA
        if (profile_enabled) decode_seconds = decode_seconds + f_poem_elapsed(started)
        call s_poem_check(code)
#else
        call s_poem_decode_batch(count, y, moments)
#endif

    end subroutine s_poem_decode_device

    subroutine s_poem_advance_device(count, active, pressure, h, y, following, ratio, root)

        integer, intent(in)        :: count
        integer(c_int), intent(in) :: active(count)
        real(wp), intent(in)       :: pressure(count), h, y(poem_dimension, count)
        real(wp), intent(out)      :: following(poem_dimension, count), ratio(count), root(count)
        integer(c_int)             :: code
        integer(int64)             :: started

        call s_poem_initialize()
#ifdef MFC_OpenACC
        !$acc wait
        source_calls = source_calls + 1_int64
        source_states = source_states + int(count, int64)
        if (profile_enabled) call system_clock(started)
        #:call GPU_HOST_DATA(use_device_addr='[active, pressure, y, following, ratio, root]')
            code = c_step_device(int(count, c_int), active, pressure, real(h, c_double), y, following, ratio, root, 0_int64)
        #:endcall GPU_HOST_DATA
        if (profile_enabled) source_seconds = source_seconds + f_poem_elapsed(started)
        call s_poem_check(code)
#else
        call s_poem_advance_active(count, active == 1_c_int, pressure, h, y, following, ratio, root)
#endif

    end subroutine s_poem_advance_device

    subroutine s_poem_decode(y, moments)

        real(wp), intent(in)  :: y(poem_dimension)
        real(wp), intent(out) :: moments(4)
        real(wp)              :: batch_y(poem_dimension, 1), batch_moments(4, 1)

        batch_y(:,1) = y
        call s_poem_decode_batch(1, batch_y, batch_moments)
        moments = batch_moments(:,1)

    end subroutine s_poem_decode

    !> One ABI call; face callers always use this uncached entry point.
    subroutine s_poem_decode_batch(count, y, moments)

        integer, intent(in)   :: count
        real(wp), intent(in)  :: y(poem_dimension, count)
        real(wp), intent(out) :: moments(4, count)
        real(c_double)        :: yc(poem_dimension, count), mc(4, count)
        integer(c_int)        :: valid(count), code
        integer(int64)        :: started

        call s_poem_initialize()
        if (count <= 0) call s_poem_abort('POEM decode requires a positive batch size.')
        yc = y
        valid = 0_c_int
        decode_calls = decode_calls + 1_int64
        decode_states = decode_states + int(count, int64)
        if (profile_enabled) call system_clock(started)
        code = c_decode(int(count, c_int), yc, mc, valid)
        if (profile_enabled) decode_seconds = decode_seconds + f_poem_elapsed(started)
        call s_poem_check(code)
        if (any(valid /= 1_c_int) .or. .not. all(ieee_is_finite(mc))) &
            & call s_poem_abort('POEM decoder returned invalid rows despite a successful status.')
        if (any(mc(1,:) <= 0._c_double)) call s_poem_abort('POEM decoder returned nonpositive bubble volume.')
        moments = mc

    end subroutine s_poem_decode_batch

    subroutine s_poem_invalidate_cells()

        if (allocated(cell_valid)) cell_valid = .false.

    end subroutine s_poem_invalidate_cells

    subroutine s_poem_size_cell_cache(count)

        integer, intent(in) :: count

        if (allocated(cell_valid)) then
            if (size(cell_valid) == count) return
            deallocate (cell_bits, cell_moments, cell_valid)
        end if
        allocate (cell_bits(poem_dimension, count), cell_moments(4, count), cell_valid(count))
        cell_valid = .false.

    end subroutine s_poem_size_cell_cache

    !> Cache only the corresponding cell's exact state input. Neither pressure nor number density enters decode; ptil is recomputed
    !! by callers. Fresh misses use power-of-two batches to bound JAX shape compilation.
    subroutine s_poem_decode_cells(count, y, moments)

        integer, intent(in)   :: count
        real(wp), intent(in)  :: y(poem_dimension, count)
        real(wp), intent(out) :: moments(4, count)
        integer               :: j, row, missing, padded, indices(count)
        integer(int64)        :: bits(poem_dimension)
        logical               :: hit
        real(wp), allocatable :: batch_y(:,:), batch_moments(:,:)

        call s_poem_initialize()
        if (count <= 0) call s_poem_abort('POEM cell decode requires a positive batch size.')
        call s_poem_size_cell_cache(count)
        missing = 0
        do j = 1, count
            bits = transfer(y(:,j), bits)
            hit = .false.
            if (cell_valid(j)) hit = all(ieee_is_finite(y(:,j))) .and. all(bits == cell_bits(:,j))
            if (hit) then
                moments(:,j) = cell_moments(:,j)
            else
                cell_valid(j) = .false.
                missing = missing + 1
                indices(missing) = j
            end if
        end do
        cell_cache_hits = cell_cache_hits + int(count - missing, int64)
        cell_cache_misses = cell_cache_misses + int(missing, int64)
        if (missing == 0) return
        padded = 1
        do while (padded < missing)
            if (padded > huge(padded)/2) call s_poem_abort('POEM decode batch size overflow.')
            padded = 2*padded
        end do
        cell_padding_states = cell_padding_states + int(padded - missing, int64)
        allocate (batch_y(poem_dimension, padded), batch_moments(4, padded))
        do row = 1, missing
            batch_y(:,row) = y(:,indices(row))
        end do
        ! Repeat a requested row, never introduce a synthetic state or clipping.
        do row = missing + 1, padded
            batch_y(:,row) = batch_y(:,1)
        end do
        call s_poem_decode_batch(padded, batch_y, batch_moments)
        do row = 1, missing
            j = indices(row)
            moments(:,j) = batch_moments(:,row)
            cell_moments(:,j) = batch_moments(:,row)
            cell_bits(:,j) = transfer(y(:,j), bits)
            cell_valid(j) = .true.
        end do

    end subroutine s_poem_decode_cells

    subroutine s_poem_advance(count, pressure, h, y, following, ratio, root)

        integer, intent(in)   :: count
        real(wp), intent(in)  :: pressure(count), h, y(poem_dimension, count)
        real(wp), intent(out) :: following(poem_dimension, count), ratio(count), root(count)
        real(c_double)        :: pc(count), yc(poem_dimension, count), nextc(poem_dimension, count), moments(4, count), rc(count), &
             & ec(count)
        integer(c_int) :: valid(count), code
        integer(int64) :: started, bits(poem_dimension)
        integer        :: j

        call s_poem_initialize()
        if (count <= 0) call s_poem_abort('POEM source requires a positive batch size.')
        call s_poem_invalidate_cells()
        pc = pressure; yc = y
        valid = 0_c_int
        source_calls = source_calls + 1_int64
        source_states = source_states + int(count, int64)
        if (profile_enabled) call system_clock(started)
        code = c_step(int(count, c_int), pc, real(h, c_double), yc, nextc, moments, valid, rc, ec)
        if (profile_enabled) source_seconds = source_seconds + f_poem_elapsed(started)
        call s_poem_check(code)
        if (any(valid /= 1_c_int) .or. .not. all(ieee_is_finite(nextc)) .or. .not. all(ieee_is_finite(moments))) &
            & call s_poem_abort('POEM source returned invalid rows despite a successful status.')
        if (any(moments(1,:) <= 0._c_double)) call s_poem_abort('POEM source returned nonpositive bubble volume.')
        following = nextc; ratio = rc; root = ec
        call s_poem_size_cell_cache(count)
        do j = 1, count
            ! Key the raw state actually decoded by the ABI, before the
            ! number*following / number conservative round trip in the host.
            cell_bits(:,j) = transfer(nextc(:,j), bits)
        end do
        cell_moments = moments
        cell_valid = .true.

    end subroutine s_poem_advance

    !> Advance only cells whose host bubble source is active. Inactive rows never enter the model. Preserve the existing all-active
    !! batch exactly; mixed batches repeat a genuine active row to bound JAX compilation shapes.
    subroutine s_poem_advance_active(count, active, pressure, h, y, following, ratio, root)

        integer, intent(in)   :: count
        logical, intent(in)   :: active(count)
        real(wp), intent(in)  :: pressure(count), h, y(poem_dimension, count)
        real(wp), intent(out) :: following(poem_dimension, count), ratio(count), root(count)
        integer               :: j, row, nactive, padded, indices(count)
        real(wp), allocatable :: batch_y(:,:), batch_pressure(:), batch_following(:,:), batch_ratio(:), batch_root(:)

        if (count <= 0) call s_poem_abort('POEM source requires a positive batch size.')
        following = y
        ratio = 0._wp
        root = 0._wp
        nactive = 0
        do j = 1, count
            if (.not. active(j)) cycle
            nactive = nactive + 1
            indices(nactive) = j
        end do
        if (nactive == 0) return
        if (nactive == count) then
            call s_poem_advance(count, pressure, h, y, following, ratio, root)
            return
        end if

        padded = 1
        do while (padded < nactive)
            if (padded > huge(padded)/2) call s_poem_abort('POEM source batch size overflow.')
            padded = 2*padded
        end do
        allocate (batch_y(poem_dimension, padded), batch_pressure(padded), batch_following(poem_dimension, padded), &
                  & batch_ratio(padded), batch_root(padded))
        do row = 1, nactive
            batch_y(:,row) = y(:,indices(row))
            batch_pressure(row) = pressure(indices(row))
        end do
        do row = nactive + 1, padded
            batch_y(:,row) = batch_y(:,1)
            batch_pressure(row) = batch_pressure(1)
        end do
        call s_poem_advance(padded, batch_pressure, h, batch_y, batch_following, batch_ratio, batch_root)
        ! Compacted row numbers are not host cell indices. Do not expose the
        ! source-populated cache to a subsequent full-cell decoder call.
        call s_poem_invalidate_cells()
        do row = 1, nactive
            j = indices(row)
            following(:,j) = batch_following(:,row)
            ratio(j) = batch_ratio(row)
            root(j) = batch_root(row)
        end do

    end subroutine s_poem_advance_active

    function f_poem_elapsed(started) result(seconds)

        integer(int64), intent(in) :: started
        integer(int64)             :: finished, rate, maximum, ticks
        real(wp)                   :: seconds

        call system_clock(finished, rate, maximum)
        ticks = finished - started
        if (finished < started) ticks = (maximum - started) + finished + 1_int64
        seconds = 0._wp
        if (rate > 0_int64) seconds = real(ticks, wp)/real(rate, wp)

    end function f_poem_elapsed

    subroutine s_poem_write_profile(path)

        character(len=*), intent(in) :: path
        integer                      :: unit

        if (.not. profile_enabled) return
        open (newunit=unit, file=path, status='replace', action='write')
        write (unit, '(A)') '# POEM process-local backend counts; wall seconds include first-use setup (JIT for JAX only).'
        write (unit, '(A,1X,I0)') 'decode_calls', decode_calls
        write (unit, '(A,1X,I0)') 'decode_states', decode_states
        write (unit, '(A,1X,I0)') 'source_calls', source_calls
        write (unit, '(A,1X,I0)') 'source_states', source_states
        write (unit, '(A,1X,I0)') 'cell_cache_hits', cell_cache_hits
        write (unit, '(A,1X,I0)') 'cell_cache_misses', cell_cache_misses
        write (unit, '(A,1X,I0)') 'cell_padding_states', cell_padding_states
        write (unit, '(A,1X,ES25.16E3)') 'decode_seconds', decode_seconds
        write (unit, '(A,1X,ES25.16E3)') 'source_seconds', source_seconds
        close (unit)

    end subroutine s_poem_write_profile

end module m_poem_backend
