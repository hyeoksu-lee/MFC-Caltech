!>
!! @file
!! @brief Contains module @ref m_bubbles_ee "m_bubbles_EE"

#:include 'macros.fpp'

!> @brief Computes ensemble-averaged (Euler--Euler) bubble source terms for radius, velocity, pressure, and mass transfer
module m_bubbles_EE

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_variables_conversion
    use m_eos
    use m_bubbles
    use m_poem_backend, only: s_poem_initialize, s_poem_decode_cells, s_poem_advance_active, poem_Aref
    use m_poem_backend, only: s_poem_decode_device, s_poem_advance_device
    use iso_c_binding, only: c_int
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none

    real(wp), allocatable, dimension(:,:,:)   :: bub_adv_src
    real(wp), allocatable, dimension(:,:,:,:) :: bub_r_src, bub_v_src, bub_p_src, bub_m_src
    $:GPU_DECLARE(create='[bub_adv_src, bub_r_src, bub_v_src, bub_p_src, bub_m_src]')

    type(scalar_field) :: divu  !< matrix for div(u)
    $:GPU_DECLARE(create='[divu]')

    integer, allocatable, dimension(:) :: rs, vs, ms, ps
    real(wp)                           :: poem_max_source_ratio = 0._wp, poem_max_root_residual = 0._wp
    $:GPU_DECLARE(create='[rs, vs, ms, ps]')
    real(wp), allocatable :: poem_gpu_y(:,:), poem_gpu_next(:,:), poem_gpu_p(:), poem_gpu_ratio(:), poem_gpu_root(:)
    real(wp), allocatable :: poem_gpu_moments(:,:)
    integer(c_int), allocatable :: poem_gpu_active(:)
    $:GPU_DECLARE(create='[poem_gpu_y, poem_gpu_next, poem_gpu_p, poem_gpu_ratio, poem_gpu_root, poem_gpu_moments, poem_gpu_active]')

contains

    !> Initialize the Euler-Euler bubble module
    impure subroutine s_initialize_bubbles_EE_module

        integer :: l

        @:ALLOCATE(rs(1:nb))
        @:ALLOCATE(vs(1:nb))
        @:ALLOCATE(ps(1:nb))
        @:ALLOCATE(ms(1:nb))

        do l = 1, nb
            rs(l) = qbmm_idx%rs(l)
            vs(l) = qbmm_idx%vs(l)
            if (.not. polytropic) then
                ps(l) = qbmm_idx%ps(l)
                ms(l) = qbmm_idx%ms(l)
            else
                ps(l) = rs(l)
                ms(l) = rs(l)
            end if
        end do

        $:GPU_UPDATE(device='[rs, vs]')
        $:GPU_UPDATE(device='[ps, ms]')

        @:ALLOCATE(divu%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))
        @:ACC_SETUP_SFs(divu)

        @:ALLOCATE(bub_adv_src(0:m, 0:n, 0:p))
        @:ALLOCATE(bub_r_src(0:m, 0:n, 0:p, 1:nb))
        @:ALLOCATE(bub_v_src(0:m, 0:n, 0:p, 1:nb))
        @:ALLOCATE(bub_p_src(0:m, 0:n, 0:p, 1:nb))
        @:ALLOCATE(bub_m_src(0:m, 0:n, 0:p, 1:nb))

        if (adap_dt .and. f_is_default(adap_dt_tol)) adap_dt_tol = dflt_adap_dt_tol
#ifdef MFC_OpenACC
        if (poem_euler) then
            @:ALLOCATE(poem_gpu_y(poem_dimension,m+1), poem_gpu_next(poem_dimension,m+1), poem_gpu_moments(4,m+1))
            @:ALLOCATE(poem_gpu_p(m+1), poem_gpu_ratio(m+1), poem_gpu_root(m+1), poem_gpu_active(m+1))
        end if
#endif

    end subroutine s_initialize_bubbles_EE_module

    !> Compute the bubble volume fraction alpha from the bubble number density
    subroutine s_comp_alpha_from_n(q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        real(wp)                                               :: nR3bar
        integer                                                :: i, j, k, l

        if (poem_euler) then
            call s_poem_initialize()
            $:GPU_PARALLEL_LOOP(private='[j]')
            do j = 0, m
                q_cons_vf(eqn_idx%alf)%sf(j, 0, 0) = (4._wp*pi/3._wp)*(poem_Aref*q_cons_vf(eqn_idx%n)%sf(j, 0, &
                          & 0) + 3._wp*q_cons_vf(eqn_idx%bub%beg)%sf(j, 0, 0))
            end do
            $:END_GPU_PARALLEL_LOOP()
            return
        end if

        $:GPU_PARALLEL_LOOP(private='[i, j, k, l, nR3bar]', collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    nR3bar = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, nb
                        nR3bar = nR3bar + weight(i)*(q_cons_vf(rs(i))%sf(j, k, l))**3._wp
                    end do
                    q_cons_vf(eqn_idx%alf)%sf(j, k, l) = (4._wp*pi*nR3bar)/(3._wp*q_cons_vf(eqn_idx%n)%sf(j, k, l)**2._wp)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_comp_alpha_from_n

    !> Frozen-pressure half source step; density, momentum, energy and number are unchanged.
    subroutine s_advance_poem_source(q_cons_vf, q_prim_vf)

        type(scalar_field), intent(inout) :: q_cons_vf(sys_size)
        type(scalar_field), intent(in)    :: q_prim_vf(sys_size)
        real(wp)                          :: y(poem_dimension, m + 1), following(poem_dimension, m + 1), pressure(m + 1), ratio(m + 1), root(m + 1), number
        logical                           :: active(m + 1)
        integer                           :: j, i

#ifdef MFC_OpenACC
        call s_advance_poem_source_gpu(q_cons_vf,q_prim_vf)
        return
#endif
        y = 0._wp
        pressure = 0._wp
        do j = 0, m
            ! Match the direct-source cutoff, including an active equality case.
            ! Re-evaluate at every half step; advection can change the local void.
            if (.not. ieee_is_finite(q_prim_vf(eqn_idx%alf)%sf(j, 0, 0))) &
                call s_mpi_abort('POEM source requires finite bubble volume fraction.')
            active(j + 1) = .not. (q_prim_vf(eqn_idx%alf)%sf(j, 0, 0) < small_alf)
            if (.not. active(j + 1)) cycle
            pressure(j + 1) = q_prim_vf(eqn_idx%E)%sf(j, 0, 0)
            do i = 1, poem_dimension
                y(i, j + 1) = q_prim_vf(eqn_idx%bub%beg + i - 1)%sf(j, 0, 0)
            end do
        end do
        call s_poem_advance_active(m + 1, active, pressure, dt/2._wp, y, following, ratio, root)
        poem_max_source_ratio = max(poem_max_source_ratio, maxval(ratio))
        poem_max_root_residual = max(poem_max_root_residual, maxval(root))
        do j = 0, m
            ! Do not rewrite inactive nY: a divide/multiply round trip need not
            ! be bitwise identity, even when the primitive source is zero.
            if (.not. active(j + 1)) cycle
            number = q_cons_vf(eqn_idx%n)%sf(j, 0, 0)
            if (number <= 0._wp) call s_mpi_abort('POEM source requires positive bubble number density.')
            do i = 1, poem_dimension
                q_cons_vf(eqn_idx%bub%beg + i - 1)%sf(j, 0, 0) = number*following(i, j + 1)
            end do
        end do

    end subroutine s_advance_poem_source

#ifdef MFC_OpenACC
    subroutine s_advance_poem_source_gpu(q_cons_vf,q_prim_vf)
        type(scalar_field), intent(inout) :: q_cons_vf(sys_size)
        type(scalar_field), intent(in) :: q_prim_vf(sys_size)
        integer :: j,i,bad
        real(wp) :: number, alf, max_ratio, max_root
        bad=0
        $:GPU_PARALLEL_LOOP(private='[j,i,alf,number]', reduction='[[bad]]', reductionOp='[max]')
        do j=0,m
            alf=q_prim_vf(eqn_idx%alf)%sf(j,0,0)
            if (.not. ieee_is_finite(alf)) bad=1
            poem_gpu_active(j+1)=0_c_int
            if (.not. (alf < small_alf)) poem_gpu_active(j+1)=1_c_int
            number=q_cons_vf(eqn_idx%n)%sf(j,0,0)
            if (poem_gpu_active(j+1)==1 .and. (.not. ieee_is_finite(number) .or. number <= 0._wp)) bad=1
            poem_gpu_p(j+1)=q_prim_vf(eqn_idx%E)%sf(j,0,0)
            $:GPU_LOOP(parallelism='[seq]')
            do i=1,poem_dimension
                poem_gpu_y(i,j+1)=q_prim_vf(eqn_idx%bub%beg+i-1)%sf(j,0,0)
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        if (bad/=0) call s_mpi_abort('Invalid POEM void fraction or active bubble number density.')
        call s_poem_advance_device(m+1,poem_gpu_active,poem_gpu_p,dt/2._wp,poem_gpu_y, &
                                  poem_gpu_next,poem_gpu_ratio,poem_gpu_root)
        max_ratio=0._wp; max_root=0._wp
        $:GPU_PARALLEL_LOOP(private='[j,i,number]', reduction='[[max_ratio,max_root]]', reductionOp='[max]')
        do j=0,m
            max_ratio=max(max_ratio,poem_gpu_ratio(j+1))
            max_root=max(max_root,poem_gpu_root(j+1))
            if (poem_gpu_active(j+1)==0) cycle
            number=q_cons_vf(eqn_idx%n)%sf(j,0,0)
            $:GPU_LOOP(parallelism='[seq]')
            do i=1,poem_dimension
                q_cons_vf(eqn_idx%bub%beg+i-1)%sf(j,0,0)=number*poem_gpu_next(i,j+1)
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        poem_max_source_ratio=max(poem_max_source_ratio,max_ratio)
        poem_max_root_residual=max(poem_max_root_residual,max_root)
    end subroutine

    subroutine s_decode_poem_cells_gpu(q_prim_vf)
        type(scalar_field), intent(in) :: q_prim_vf(sys_size)
        integer :: j,i
        $:GPU_PARALLEL_LOOP(private='[j,i]', collapse=2)
        do j=0,m
            do i=1,poem_dimension
                poem_gpu_y(i,j+1)=q_prim_vf(eqn_idx%bub%beg+i-1)%sf(j,0,0)
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
        call s_poem_decode_device(m+1,poem_gpu_y,poem_gpu_moments)
    end subroutine
#endif

    !> Compact diagnostics at full-step endpoints for the one-dimensional CPU study.
    subroutine s_write_poem_diagnostics(q_cons_vf, q_prim_vf, step, time)

        type(scalar_field), intent(in) :: q_cons_vf(sys_size), q_prim_vf(sys_size)
        integer, intent(in)            :: step
        real(wp), intent(in)           :: time
        integer                        :: unit, j, i, q
        real(wp)                       :: moments(4), radius, velocity, number
        real(wp), allocatable          :: y(:,:), decoded(:,:)
        character(len=1024)            :: path

        if (poem_euler) then
            if (n /= 0 .or. p /= 0) call s_mpi_abort('POEM diagnostics require the supported one-dimensional domain.')
            allocate(y(poem_dimension, m + 1), decoded(4, m + 1))
#ifdef MFC_OpenACC
            call s_decode_poem_cells_gpu(q_prim_vf)
            $:GPU_UPDATE(host='[poem_gpu_moments]')
            decoded=poem_gpu_moments
#else
            do j = 0, m
                do i = 1, poem_dimension
                    y(i, j + 1) = q_prim_vf(eqn_idx%bub%beg + i - 1)%sf(j, 0, 0)
                end do
            end do
            call s_poem_decode_cells(m + 1, y, decoded)
#endif
        end if
#ifdef MFC_OpenACC
        ! Explicit output-only host transfer; simulation kernels stay on device.
        do i=1,sys_size
            $:GPU_UPDATE(host='[q_cons_vf(i)%sf, q_prim_vf(i)%sf]')
        end do
#endif
        path = trim(case_dir) // '/D/poem_moments.dat'
        if (step == t_step_start) then
            open (newunit=unit, file=trim(path), status='replace', action='write')
            write (unit, '(A)') '# step time cell x rho u p alpha n A K W Q max_source_energy_ratio max_root_residual'
        else
            open (newunit=unit, file=trim(path), status='old', position='append', action='write')
        end if
        do j = 0, m
            if (poem_euler) then
                moments = decoded(:, j + 1)
            else
                moments = 0._wp
                do q = 1, nb
                    radius = q_prim_vf(rs(q))%sf(j, 0, 0)
                    velocity = q_prim_vf(vs(q))%sf(j, 0, 0)
                    moments(1) = moments(1) + weight(q)*radius**3
                    moments(2) = moments(2) + weight(q)*radius**3*velocity**2
                    moments(3) = moments(3) + weight(q)*radius**3*f_cpbw_KM(R0(q), radius, velocity, 0._wp)
                    moments(4) = moments(4) + weight(q)*radius**2*velocity
                end do
            end if
            number = q_cons_vf(eqn_idx%n)%sf(j, 0, 0)
            write (unit, '(I9,1X,ES25.16E3,1X,I9,12(1X,ES25.16E3))') step, time, j, x_cc(j), q_cons_vf(1)%sf(j, 0, 0), &
                   & q_prim_vf(eqn_idx%mom%beg)%sf(j, 0, 0), q_prim_vf(eqn_idx%E)%sf(j, 0, 0), q_prim_vf(eqn_idx%alf)%sf(j, 0, &
                   & 0), number, moments, poem_max_source_ratio, poem_max_root_residual
        end do
        close (unit)

    end subroutine s_write_poem_diagnostics

    !> Compute the right-hand side for Euler-Euler bubble transport
    subroutine s_compute_bubbles_EE_rhs(idir, q_prim_vf, divu_in)

        integer, intent(in)                                 :: idir
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), intent(inout)                   :: divu_in  !< matrix for div(u)
        integer                                             :: j, k, l

        if (idir == 1) then
            if (.not. qbmm) then
                $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            divu_in%sf(j, k, l) = 0._wp
                            divu_in%sf(j, k, l) = 5.e-1_wp/dx(j)*(q_prim_vf(eqn_idx%cont%end + idir)%sf(j + 1, k, &
                                       & l) - q_prim_vf(eqn_idx%cont%end + idir)%sf(j - 1, k, l))
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        else if (idir == 2) then
            $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        divu_in%sf(j, k, l) = divu_in%sf(j, k, l) + 5.e-1_wp/dy(k)*(q_prim_vf(eqn_idx%cont%end + idir)%sf(j, &
                                   & k + 1, l) - q_prim_vf(eqn_idx%cont%end + idir)%sf(j, k - 1, l))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        else if (idir == 3) then
            $:GPU_PARALLEL_LOOP(private='[j, k, l]', collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        divu_in%sf(j, k, l) = divu_in%sf(j, k, l) + 5.e-1_wp/dz(l)*(q_prim_vf(eqn_idx%cont%end + idir)%sf(j, k, &
                                   & l + 1) - q_prim_vf(eqn_idx%cont%end + idir)%sf(j, k, l - 1))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_compute_bubbles_EE_rhs

    !> Compute the Euler-Euler bubble source terms
    impure subroutine s_compute_bubble_EE_source(q_cons_vf, q_prim_vf, rhs_vf, divu_in)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: rhs_vf
        type(scalar_field), intent(in) :: divu_in  !< matrix for div(u)
        real(wp) :: rddot
        real(wp) :: pb_local, mv_local, vflux, pbdot
        real(wp) :: n_tait, B_tait, qv_bub
        real(wp) :: chi_vw_l, k_mw_l, rho_mw_l     !< Per-thread bubble-wall scratch (avoid module-scalar race)

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3) :: Rtmp, Vtmp
            real(wp), dimension(3) :: myalpha, myalpha_rho
        #:else
            real(wp), dimension(nb)         :: Rtmp, Vtmp
            real(wp), dimension(num_fluids) :: myalpha, myalpha_rho
        #:endif
        real(wp) :: myR, myV, alf, myP, myRho, R2Vav, R3
        real(wp) :: nbub                            !< Bubble number density
        integer  :: i, j, k, l, q, ii               !< Loop variables
        integer  :: adap_dt_stop_sum, adap_dt_stop  !< Fail-safe exit if max iteration count reached
        integer  :: dmBub_id                        !< Dummy variables for unified subgrid bubble subroutines
        real(wp) :: dmMass_v, dmMass_n, dmBeta_c, dmBeta_t, dmCson

        if (poem_euler) then
            call s_advance_poem_source(q_cons_vf, q_prim_vf)
            return
        end if

        $:GPU_PARALLEL_LOOP(private='[j, k, l, q]', collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    bub_adv_src(j, k, l) = 0._wp

                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, nb
                        bub_r_src(j, k, l, q) = 0._wp
                        bub_v_src(j, k, l, q) = 0._wp
                        bub_p_src(j, k, l, q) = 0._wp
                        bub_m_src(j, k, l, q) = 0._wp
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        adap_dt_stop_sum = 0
        $:GPU_PARALLEL_LOOP(private='[j, k, l, Rtmp, Vtmp, myalpha_rho, myalpha, myR, myV, alf, myP, myRho, R2Vav, R3, nbub, &
                            & pb_local, mv_local, vflux, pbdot, rddot, n_tait, B_tait, qv_bub, adap_dt_stop, chi_vw_l, k_mw_l, &
                            & rho_mw_l]', collapse=3, copy='[adap_dt_stop_sum]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    if (adv_n) then
                        nbub = q_prim_vf(eqn_idx%n)%sf(j, k, l)
                    else
                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, nb
                            Rtmp(q) = q_prim_vf(rs(q))%sf(j, k, l)
                            Vtmp(q) = q_prim_vf(vs(q))%sf(j, k, l)
                        end do

                        R3 = 0._wp

                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, nb
                            R3 = R3 + weight(q)*Rtmp(q)**3._wp
                        end do

                        nbub = (3._wp/(4._wp*pi))*q_prim_vf(eqn_idx%alf)%sf(j, k, l)/R3
                    end if

                    if (.not. adap_dt) then
                        R2Vav = 0._wp

                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, nb
                            R2Vav = R2Vav + weight(q)*Rtmp(q)**2._wp*Vtmp(q)
                        end do

                        bub_adv_src(j, k, l) = 4._wp*pi*nbub*R2Vav
                    end if

                    $:GPU_LOOP(parallelism='[seq]')
                    do q = 1, nb
                        $:GPU_LOOP(parallelism='[seq]')
                        do ii = 1, num_fluids
                            myalpha_rho(ii) = q_cons_vf(ii)%sf(j, k, l)
                            myalpha(ii) = q_cons_vf(eqn_idx%adv%beg + ii - 1)%sf(j, k, l)
                        end do

                        call s_compute_mixture_coefficients(myalpha_rho, myalpha, myRho, n_tait, B_tait, qv_bub)
                        B_tait = B_tait/pi_fac

                        ! n_tait still holds the mixture gamma here, which the reference pressure needs
                        B_tait = f_isentrope_pressure(B_tait, n_tait)
                        n_tait = f_isentrope_exponent(n_tait)

                        myP = q_prim_vf(eqn_idx%E)%sf(j, k, l)
                        alf = q_prim_vf(eqn_idx%alf)%sf(j, k, l)
                        myR = q_prim_vf(rs(q))%sf(j, k, l)
                        myV = q_prim_vf(vs(q))%sf(j, k, l)

                        if (alf < small_alf) then
                            bub_adv_src(j, k, l) = 0._wp
                            bub_r_src(j, k, l, q) = 0._wp
                            bub_v_src(j, k, l, q) = 0._wp
                            if (.not. polytropic) then
                                bub_p_src(j, k, l, q) = 0._wp
                                bub_m_src(j, k, l, q) = 0._wp
                            end if
                        else
                            if (.not. polytropic) then
                                pb_local = q_prim_vf(ps(q))%sf(j, k, l)
                                mv_local = q_prim_vf(ms(q))%sf(j, k, l)
                                call s_bwproperty(pb_local, q, chi_vw_l, k_mw_l, rho_mw_l)
                                call s_vflux(myR, myV, pb_local, mv_local, q, vflux, fchi_vw=chi_vw_l, frho_mw=rho_mw_l)
                                pbdot = f_bpres_dot(vflux, myR, myV, pb_local, mv_local, q, fk_mw=k_mw_l)
                                bub_p_src(j, k, l, q) = nbub*pbdot
                                bub_m_src(j, k, l, q) = nbub*vflux*4._wp*pi*(myR**2._wp)
                            else
                                pb_local = 0._wp; mv_local = 0._wp; vflux = 0._wp; pbdot = 0._wp
                            end if

                            adap_dt_stop = 0

                            ! Adaptive time stepping
                            if (adap_dt) then
                                ! Matched POEM benchmark only; flow density/void fraction remain untouched.
                                if (poem_fixed_density) then
                                    myRho = 1._wp
                                    alf = 0._wp
                                end if
                                adap_dt_stop = f_advance_step(myRho, myP, myR, myV, R0(q), pb_local, pbdot, alf, n_tait, B_tait, &
                                                              & bub_adv_src(j, k, l), divu_in%sf(j, k, l), dmBub_id, dmMass_v, &
                                                              & dmMass_n, dmBeta_c, dmBeta_t, dmCson)

                                q_cons_vf(rs(q))%sf(j, k, l) = nbub*myR
                                q_cons_vf(vs(q))%sf(j, k, l) = nbub*myV
                            else
                                rddot = f_rddot(myRho, myP, myR, myV, R0(q), pb_local, pbdot, alf, n_tait, B_tait, bub_adv_src(j, &
                                                & k, l), divu_in%sf(j, k, l), dmCson)
                                bub_v_src(j, k, l, q) = nbub*rddot
                                bub_r_src(j, k, l, q) = q_cons_vf(vs(q))%sf(j, k, l)
                            end if

                            $:GPU_ATOMIC(atomic='update')
                            adap_dt_stop_sum = adap_dt_stop_sum + adap_dt_stop
                        end if
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        if (adap_dt .and. adap_dt_stop_sum > 0) call s_mpi_abort("Adaptive time stepping failed to converge.")

        if (.not. adap_dt) then
            $:GPU_PARALLEL_LOOP(private='[i, k, l, q]', collapse=3)
            do l = 0, p
                do q = 0, n
                    do i = 0, m
                        rhs_vf(eqn_idx%alf)%sf(i, q, l) = rhs_vf(eqn_idx%alf)%sf(i, q, l) + bub_adv_src(i, q, l)
                        if (num_fluids > 1) rhs_vf(eqn_idx%adv%beg)%sf(i, q, l) = rhs_vf(eqn_idx%adv%beg)%sf(i, q, &
                            & l) - bub_adv_src(i, q, l)
                        $:GPU_LOOP(parallelism='[seq]')
                        do k = 1, nb
                            rhs_vf(rs(k))%sf(i, q, l) = rhs_vf(rs(k))%sf(i, q, l) + bub_r_src(i, q, l, k)
                            rhs_vf(vs(k))%sf(i, q, l) = rhs_vf(vs(k))%sf(i, q, l) + bub_v_src(i, q, l, k)
                            if (polytropic .neqv. .true.) then
                                rhs_vf(ps(k))%sf(i, q, l) = rhs_vf(ps(k))%sf(i, q, l) + bub_p_src(i, q, l, k)
                                rhs_vf(ms(k))%sf(i, q, l) = rhs_vf(ms(k))%sf(i, q, l) + bub_m_src(i, q, l, k)
                            end if
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_compute_bubble_EE_source

    !> Evaluate the Euler-Euler pressure correction at cell centers and store it in ptil. This mirrors the modified mixture pressure
    !! that the HLLC solver applies to the reconstructed face states in m_riemann_solver_hllc.fpp, following Ando (2010), Eq. (2.5).
    !! It is a diagnostic: ptil is read only by the probe output path and never feeds back into the solution.
    impure subroutine s_compute_ptilde(q_cons_vf, q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_cons_vf, q_prim_vf
        real(wp)                                            :: R3bar, R3V2bar, PbwR3bar
        real(wp)                                            :: myRho, myP, alf, myR, myV, myPb
        real(wp), allocatable                              :: poem_y(:,:), poem_moments(:,:)
        integer                                             :: j, k, l, q, ii

#ifdef MFC_OpenACC
        if (poem_euler) call s_decode_poem_cells_gpu(q_prim_vf)
#endif
#if !defined(MFC_OpenACC) && !defined(MFC_OpenMP)
        if (poem_euler) then
            if (n /= 0 .or. p /= 0) call s_mpi_abort('POEM cell decoding requires the supported one-dimensional domain.')
            allocate(poem_y(poem_dimension, m + 1), poem_moments(4, m + 1))
            do j = 0, m
                do q = 1, poem_dimension
                    poem_y(q, j + 1) = q_prim_vf(eqn_idx%bub%beg + q - 1)%sf(j, 0, 0)
                end do
            end do
            call s_poem_decode_cells(m + 1, poem_y, poem_moments)
        end if
#endif
        $:GPU_PARALLEL_LOOP(private='[j, k, l, q, ii, R3bar, R3V2bar, PbwR3bar, myRho, myP, alf, myR, myV, myPb]', collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    alf = q_prim_vf(eqn_idx%alf)%sf(j, k, l)
                    myP = q_prim_vf(eqn_idx%E)%sf(j, k, l)

                    myRho = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do ii = 1, num_fluids
                        myRho = myRho + q_cons_vf(ii)%sf(j, k, l)
                    end do

                    if (poem_euler) then
#ifdef MFC_OpenACC
                        R3bar = poem_gpu_moments(1,j+1)
                        R3V2bar = poem_gpu_moments(2,j+1)
                        PbwR3bar = poem_gpu_moments(3,j+1)
#endif
#if !defined(MFC_OpenACC) && !defined(MFC_OpenMP)
                        R3bar = poem_moments(1, j + 1)
                        R3V2bar = poem_moments(2, j + 1)
                        PbwR3bar = poem_moments(3, j + 1)
#endif
                    else if (qbmm) then
                        R3bar = mom_sp(1)%sf(j, k, l)
                        R3V2bar = mom_sp(3)%sf(j, k, l)
                        PbwR3bar = mom_sp(4)%sf(j, k, l)
                    else
                        R3bar = 0._wp
                        R3V2bar = 0._wp
                        PbwR3bar = 0._wp

                        $:GPU_LOOP(parallelism='[seq]')
                        do q = 1, nb
                            myR = q_prim_vf(rs(q))%sf(j, k, l)
                            myV = q_prim_vf(vs(q))%sf(j, k, l)
                            if (polytropic) then
                                myPb = 0._wp
                            else
                                myPb = q_prim_vf(ps(q))%sf(j, k, l)
                            end if

                            R3bar = R3bar + weight(q)*(myR**3._wp)
                            R3V2bar = R3V2bar + weight(q)*(myR**3._wp)*(myV**2._wp)
                            PbwR3bar = PbwR3bar + weight(q)*f_cpbw_KM(R0(q), myR, myV, myPb)*(myR**3._wp)
                        end do
                    end if

                    ! Same degenerate-state guard the HLLC solver uses
                    if (alf < small_alf .or. R3bar < small_alf) then
                        ptil(j, k, l) = alf*myP
                    else
                        ptil(j, k, l) = alf*(myP - PbwR3bar/R3bar - myRho*R3V2bar/R3bar)
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_ptilde

end module m_bubbles_EE
