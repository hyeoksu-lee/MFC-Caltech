!>
!!@file
!!@brief Contains module m_checker_common

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Shared input validation checks for grid dimensions and AMD GPU compiler limits
module m_checker_common

    use m_global_parameters
    use m_mpi_proxy
    use m_helper_basic
    use m_helper

    implicit none

    private; public :: s_check_inputs_common

contains

    !> Checks compatibility of parameters in the input file. Used by all three stages
    impure subroutine s_check_inputs_common(check_total_cells, n_global)

        logical, intent(in)         :: check_total_cells
        integer(kind=8), intent(in) :: n_global

        if (check_total_cells) call s_check_total_cells(n_global)
        if (poem_euler .or. poem_diagnostics) then
            ! lint: runtime-check MPI world size is determined by the executable launch.
            @:PROHIBIT(num_procs /= 1, "POEM research integration requires one MPI rank")
#if defined(MFC_OpenMP) || defined(MFC_SINGLE_PRECISION) || defined(MFC_MIXED_PRECISION)
            call s_mpi_abort("POEM requires CPU or NVIDIA OpenACC double precision")
#endif
        end if
#ifdef MFC_POST_PROCESS
        if (poem_euler) then
            call s_mpi_abort("POEM post_process radius-bin output is unsupported; use poem_diagnostics")
        end if
#endif
        #:if USING_AMD
            call s_check_amd
        #:endif

    end subroutine s_check_inputs_common

    !> Verify that the total number of grid cells meets the minimum required by the number of dimensions and MPI ranks.
    impure subroutine s_check_total_cells(n_global)

        character(len=18)           :: numStr  !< for int to string conversion
        integer(kind=8)             :: min_cells
        integer(kind=8), intent(in) :: n_global

        min_cells = int(2, kind=8)**int(min(1, m) + min(1, n) + min(1, p), kind=8)*int(num_procs, kind=8)
        call s_int_to_str(2**(min(1, m) + min(1, n) + min(1, p))*num_procs, numStr)

        @:PROHIBIT(n_global < min_cells, &
                   & "Total number of cells must be at least (2^[number of dimensions])*num_procs, " // "which is currently " &
                   & // trim(numStr))

    end subroutine s_check_total_cells

    !> Check that simulation parameters stay within AMD GPU compiler limits when case optimization is disabled.
    impure subroutine s_check_amd

        #:if not MFC_CASE_OPTIMIZATION
            @:PROHIBIT(num_fluids > 3, "num_fluids <= 3 for AMDFLang when Case optimization is off")
            @:PROHIBIT((bubbles_euler .or. bubbles_lagrange) .and. nb > 3, "nb <= 3 for AMDFLang when Case optimization is off")
            ! HLLC's star states and the CBC L vectors have no other bound.
            @:PROHIBIT(sys_size > ${AMD_SYS_SIZE_MAX}$, &
                       & "sys_size <= ${AMD_SYS_SIZE_MAX}$ for AMDFLang when Case optimization is off")
        #:endif

    end subroutine s_check_amd

end module m_checker_common
