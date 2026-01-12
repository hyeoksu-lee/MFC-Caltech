!>
!! @file m_time_steppers.f90
!! @brief Contains module m_time_steppers

#:include 'macros.fpp'

!> @brief The following module features a variety of time-stepping schemes.
!!              Currently, it includes the following Runge-Kutta (RK) algorithms:
!!                   1) 1st Order TVD RK
!!                   2) 2nd Order TVD RK
!!                   3) 3rd Order TVD RK
!!              where TVD designates a total-variation-diminishing time-stepper.
module m_time_steppers

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_rhs                  !< Right-hane-side (RHS) evaluation procedures

    use m_pressure_relaxation  !< Pressure relaxation procedures

    use m_data_output          !< Run-time info & solution data output procedures

    use m_bubbles_EE           !< Ensemble-averaged bubble dynamics routines

    use m_bubbles_EL           !< Lagrange bubble dynamics routines

    use m_ibm

    use m_hyperelastic

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_boundary_common

    use m_helper

    use m_sim_helpers

    use m_fftw

    use m_nvtx

    use m_thermochem, only: num_species

    use m_body_forces

    implicit none

    type(vector_field), allocatable, dimension(:) :: q_cons_ts !<
    !! Cell-average conservative variables at each time-stage (TS)

    type(vector_field), allocatable, dimension(:) :: q_cons_pts !<
    !! Cell-average conservative variables at each pseudo-time-stage (PTS)

    type(scalar_field), allocatable, dimension(:) :: q_prim_vf !<
    !! Cell-average primitive variables at the current time-stage

    type(scalar_field), allocatable, dimension(:) :: rhs_vf !<
    !! Cell-average RHS variables at the current time-stage

    type(integer_field), allocatable, dimension(:, :) :: bc_type !<
    !! Boundary condition identifiers

    type(vector_field), allocatable, dimension(:) :: q_prim_ts1, q_prim_ts2 !<
    !! Cell-average primitive variables at consecutive TIMESTEPS

    real(wp), allocatable, dimension(:, :, :, :, :) :: rhs_pb

    type(scalar_field) :: q_T_sf !<
    !! Cell-average temperature variables at the current time-stage

    real(wp), allocatable, dimension(:, :, :, :, :) :: rhs_mv

    real(wp), allocatable, dimension(:, :, :) :: max_dt

    integer, private :: num_ts !<
    !! Number of time stages in the time-stepping scheme

    real(wp), allocatable, dimension(:, :, :, :, :) :: fjacobian

    integer :: stor !< storage index
    real(wp), allocatable, dimension(:, :) :: rk_coef
    integer, private :: num_probe_ts

    $:GPU_DECLARE(create='[q_cons_ts,q_cons_pts,q_prim_vf,q_T_sf,rhs_vf,q_prim_ts,rhs_mv,rhs_pb,max_dt,rk_coef]')

#if defined(__NVCOMPILER_GPU_UNIFIED_MEM)
    real(wp), allocatable, dimension(:, :, :, :), pinned, target :: q_cons_ts_pool_host, q_cons_pts_pool_host
#elif defined(FRONTIER_UNIFIED)
    real(wp), pointer, contiguous, dimension(:, :, :, :) :: q_cons_ts_pool_host, q_cons_ts_pool_device
    real(wp), pointer, contiguous, dimension(:, :, :, :) :: q_cons_pts_pool_host, q_cons_pts_pool_device
    integer(kind=8) :: pool_dims(4), pool_starts(4)
    integer(kind=8) :: pool_size
    type(c_ptr) :: cptr_host, cptr_device
#endif

    real(wp), allocatable, dimension(:, :, :) :: rho_avg_x, vel_avg_rms_x, H_avg_x, gamma_avg_x
    real(wp), allocatable, dimension(:, :, :, :) :: vel_avg_x
    real(wp), allocatable, dimension(:, :, :) :: rho_avg_y, vel_avg_rms_y, H_avg_y, gamma_avg_y
    real(wp), allocatable, dimension(:, :, :, :) :: vel_avg_y

contains

    !> The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    impure subroutine s_initialize_time_steppers_module
#ifdef FRONTIER_UNIFIED
        use hipfort
        use hipfort_hipmalloc
        use hipfort_check
#if defined(MFC_OpenACC)
        use openacc
#endif
#endif
        integer :: i, j !< Generic loop iterators

        ! Setting number of time-stages for selected time-stepping scheme
        if (time_stepper == 1) then
            num_ts = 1
        elseif (any(time_stepper == (/2, 3, 4/))) then
            num_ts = 2
        end if

        if (probe_wrt) then
            num_probe_ts = 2
        end if

        ! Allocating the cell-average conservative variables
        @:ALLOCATE(q_cons_ts(1:num_ts))
        @:PREFER_GPU(q_cons_ts)
        @:ALLOCATE(q_cons_pts(1:num_ts))
        @:PREFER_GPU(q_cons_pts)

        do i = 1, num_ts
            @:ALLOCATE(q_cons_ts(i)%vf(1:sys_size))
            @:PREFER_GPU(q_cons_ts(i)%vf)
            @:ALLOCATE(q_cons_pts(i)%vf(1:sys_size))
            @:PREFER_GPU(q_cons_pts(i)%vf)
        end do

#if defined(__NVCOMPILER_GPU_UNIFIED_MEM)
        if (num_ts == 2 .and. nv_uvm_out_of_core) then
            ! host allocation for q_cons_ts(2)%vf(j)%sf for all j
            allocate (q_cons_ts_pool_host(idwbuff(1)%beg:idwbuff(1)%end, &
                                          idwbuff(2)%beg:idwbuff(2)%end, &
                                          idwbuff(3)%beg:idwbuff(3)%end, &
                                          1:sys_size))
            ! host allocation for q_cons_pts(2)%vf(j)%sf for all j
            allocate (q_cons_pts_pool_host(idwbuff(1)%beg:idwbuff(1)%end, &
                                           idwbuff(2)%beg:idwbuff(2)%end, &
                                           idwbuff(3)%beg:idwbuff(3)%end, &
                                           1:sys_size))
        end if

        do j = 1, sys_size
            ! q_cons_ts(1) lives on the device
            @:ALLOCATE(q_cons_ts(1)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end))
            @:PREFER_GPU(q_cons_ts(1)%vf(j)%sf)
            ! q_cons_pts(1) lives on the device
            @:ALLOCATE(q_cons_pts(1)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end))
            @:PREFER_GPU(q_cons_pts(1)%vf(j)%sf)
            if (num_ts == 2) then
                if (nv_uvm_out_of_core) then
                    ! q_cons_ts(2) lives on the host
                    q_cons_ts(2)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                                          idwbuff(2)%beg:idwbuff(2)%end, &
                                          idwbuff(3)%beg:idwbuff(3)%end) => q_cons_ts_pool_host(:, :, :, j)
                    ! q_cons_pts(2) lives on the host
                    q_cons_pts(2)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                                           idwbuff(2)%beg:idwbuff(2)%end, &
                                           idwbuff(3)%beg:idwbuff(3)%end) => q_cons_pts_pool_host(:, :, :, j)
                else
                    @:ALLOCATE(q_cons_ts(2)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                        idwbuff(2)%beg:idwbuff(2)%end, &
                        idwbuff(3)%beg:idwbuff(3)%end))
                    @:PREFER_GPU(q_cons_ts(2)%vf(j)%sf)
                    @:ALLOCATE(q_cons_pts(2)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                        idwbuff(2)%beg:idwbuff(2)%end, &
                        idwbuff(3)%beg:idwbuff(3)%end))
                    @:PREFER_GPU(q_cons_pts(2)%vf(j)%sf)
                end if
            end if
        end do

        do i = 1, num_ts
            @:ACC_SETUP_VFs(q_cons_ts(i))
            @:ACC_SETUP_VFs(q_cons_pts(i))
        end do
#elif defined(FRONTIER_UNIFIED)
        ! Allocate to memory regions using hip calls
        ! that we will attach pointers to
        do i = 1, 3
            pool_dims(i) = idwbuff(i)%end - idwbuff(i)%beg + 1
            pool_starts(i) = idwbuff(i)%beg
        end do
        pool_dims(4) = sys_size
        pool_starts(4) = 1
#ifdef MFC_MIXED_PRECISION
        pool_size = 1_8*(idwbuff(1)%end - idwbuff(1)%beg + 1)*(idwbuff(2)%end - idwbuff(2)%beg + 1)*(idwbuff(3)%end - idwbuff(3)%beg + 1)*sys_size
        call hipCheck(hipMalloc_(cptr_device, pool_size*2_8))
        call c_f_pointer(cptr_device, q_cons_ts_pool_device, shape=pool_dims)
        q_cons_ts_pool_device(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:) => q_cons_ts_pool_device

        call hipCheck(hipMallocManaged_(cptr_host, pool_size*2_8, hipMemAttachGlobal))
        call c_f_pointer(cptr_host, q_cons_ts_pool_host, shape=pool_dims)
        q_cons_ts_pool_host(idwbuff(1)%beg:, idwbuff(2)%beg:, idwbuff(3)%beg:, 1:) => q_cons_ts_pool_host
#else
        ! Doing hipMalloc then mapping should be most performant
        call hipCheck(hipMalloc(q_cons_ts_pool_device, dims8=pool_dims, lbounds8=pool_starts))
        ! Without this map CCE will still create a device copy, because it's silly like that
#if defined(MFC_OpenACC)
        call acc_map_data(q_cons_ts_pool_device, c_loc(q_cons_ts_pool_device), c_sizeof(q_cons_ts_pool_device))
#endif
        ! CCE see it can access this and will leave it on the host. It will stay on the host so long as HSA_XNACK=1
        ! NOTE: WE CANNOT DO ATOMICS INTO THIS MEMORY. We have to change a property to use atomics here
        ! Otherwise leaving this as fine-grained will actually help performance since it can't be cached in GPU L2
        if (num_ts == 2) then
            call hipCheck(hipMallocManaged(q_cons_ts_pool_host, dims8=pool_dims, lbounds8=pool_starts, flags=hipMemAttachGlobal))
#if defined(MFC_OpenMP)
            call hipCheck(hipMemAdvise(c_loc(q_cons_ts_pool_host), c_sizeof(q_cons_ts_pool_host), hipMemAdviseSetPreferredLocation, -1))
#endif
        end if
#endif

        do j = 1, sys_size
            ! q_cons_ts(1) lives on the device
            q_cons_ts(1)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                                  idwbuff(2)%beg:idwbuff(2)%end, &
                                  idwbuff(3)%beg:idwbuff(3)%end) => q_cons_ts_pool_device(:, :, :, j)
            ! q_cons_pts(1) lives on the device
            q_cons_pts(1)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                                   idwbuff(2)%beg:idwbuff(2)%end, &
                                   idwbuff(3)%beg:idwbuff(3)%end) => q_cons_pts_pool_device(:, :, :, j)
            if (num_ts == 2) then
                ! q_cons_ts(2) lives on the host
                q_cons_ts(2)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                                      idwbuff(2)%beg:idwbuff(2)%end, &
                                      idwbuff(3)%beg:idwbuff(3)%end) => q_cons_ts_pool_host(:, :, :, j)
                ! q_cons_pts(2) lives on the host
                q_cons_pts(2)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                                       idwbuff(2)%beg:idwbuff(2)%end, &
                                       idwbuff(3)%beg:idwbuff(3)%end) => q_cons_pts_pool_host(:, :, :, j)
            end if
        end do

        do i = 1, num_ts
            @:ACC_SETUP_VFs(q_cons_ts(i))
            @:ACC_SETUP_VFs(q_cons_pts(i))
            do j = 1, sys_size
                $:GPU_UPDATE(device='[q_cons_ts(i)%vf(j)]')
                $:GPU_UPDATE(device='[q_cons_pts(i)%vf(j)]')
            end do
        end do
#else
        do i = 1, num_ts
            do j = 1, sys_size
                @:ALLOCATE(q_cons_ts(i)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
                @:ALLOCATE(q_cons_pts(i)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
            end do
            @:ACC_SETUP_VFs(q_cons_ts(i))
            @:ACC_SETUP_VFs(q_cons_pts(i))
        end do
#endif

        ! Allocating the cell-average primitive ts variables
        if (probe_wrt) then
            @:ALLOCATE(q_prim_ts1(1:num_probe_ts))

            do i = 1, num_probe_ts
                @:ALLOCATE(q_prim_ts1(i)%vf(1:sys_size))
            end do

            do i = 1, num_probe_ts
                do j = 1, sys_size
                    @:ALLOCATE(q_prim_ts1(i)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                        idwbuff(2)%beg:idwbuff(2)%end, &
                        idwbuff(3)%beg:idwbuff(3)%end))
                end do
                @:ACC_SETUP_VFs(q_prim_ts1(i))
            end do

            @:ALLOCATE(q_prim_ts2(1:num_probe_ts))

            do i = 1, num_probe_ts
                @:ALLOCATE(q_prim_ts2(i)%vf(1:sys_size))
            end do

            do i = 1, num_probe_ts
                do j = 1, sys_size
                    @:ALLOCATE(q_prim_ts2(i)%vf(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                        idwbuff(2)%beg:idwbuff(2)%end, &
                        idwbuff(3)%beg:idwbuff(3)%end))
                end do
                @:ACC_SETUP_VFs(q_prim_ts2(i))
            end do
        end if

        ! Allocating the cell-average primitive variables
        @:ALLOCATE(q_prim_vf(1:sys_size))

        if (.not. igr) then
            do i = 1, adv_idx%end
                @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(q_prim_vf(i))
            end do

            if (bubbles_euler) then
                do i = bub_idx%beg, bub_idx%end
                    @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                        idwbuff(2)%beg:idwbuff(2)%end, &
                        idwbuff(3)%beg:idwbuff(3)%end))
                    @:ACC_SETUP_SFs(q_prim_vf(i))
                end do
                if (adv_n) then
                    @:ALLOCATE(q_prim_vf(n_idx)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                        idwbuff(2)%beg:idwbuff(2)%end, &
                        idwbuff(3)%beg:idwbuff(3)%end))
                    @:ACC_SETUP_SFs(q_prim_vf(n_idx))
                end if
            end if

            if (mhd) then
                do i = B_idx%beg, B_idx%end
                    @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                        idwbuff(2)%beg:idwbuff(2)%end, &
                        idwbuff(3)%beg:idwbuff(3)%end))
                    @:ACC_SETUP_SFs(q_prim_vf(i))
                end do
            end if

            if (elasticity) then
                do i = stress_idx%beg, stress_idx%end
                    @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                        idwbuff(2)%beg:idwbuff(2)%end, &
                        idwbuff(3)%beg:idwbuff(3)%end))
                    @:ACC_SETUP_SFs(q_prim_vf(i))
                end do
            end if

            if (hyperelasticity) then
                do i = xibeg, xiend + 1
                    @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                        idwbuff(2)%beg:idwbuff(2)%end, &
                        idwbuff(3)%beg:idwbuff(3)%end))
                    @:ACC_SETUP_SFs(q_prim_vf(i))
                end do
            end if

            if (cont_damage) then
                @:ALLOCATE(q_prim_vf(damage_idx)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(q_prim_vf(damage_idx))
            end if

            if (model_eqns == 3) then
                do i = internalEnergies_idx%beg, internalEnergies_idx%end
                    @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                        idwbuff(2)%beg:idwbuff(2)%end, &
                        idwbuff(3)%beg:idwbuff(3)%end))
                    @:ACC_SETUP_SFs(q_prim_vf(i))
                end do
            end if

            if (surface_tension) then
                @:ALLOCATE(q_prim_vf(c_idx)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(q_prim_vf(c_idx))
            end if

            if (chemistry) then
                do i = chemxb, chemxe
                    @:ALLOCATE(q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                        idwbuff(2)%beg:idwbuff(2)%end, &
                        idwbuff(3)%beg:idwbuff(3)%end))
                    @:ACC_SETUP_SFs(q_prim_vf(i))
                end do

                @:ALLOCATE(q_T_sf%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                    idwbuff(2)%beg:idwbuff(2)%end, &
                    idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(q_T_sf)
            end if
        end if

        @:ALLOCATE(pb_ts(1:2))
        !Initialize bubble variables pb and mv at all quadrature nodes for all R0 bins
        if (qbmm .and. (.not. polytropic)) then
            @:ALLOCATE(pb_ts(1)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(pb_ts(1))

            @:ALLOCATE(pb_ts(2)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(pb_ts(2))

            @:ALLOCATE(rhs_pb(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end, 1:nnode, 1:nb))
        else if (qbmm .and. polytropic) then
            @:ALLOCATE(pb_ts(1)%sf(idwbuff(1)%beg:idwbuff(1)%beg + 1, &
                idwbuff(2)%beg:idwbuff(2)%beg + 1, &
                idwbuff(3)%beg:idwbuff(3)%beg + 1, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(pb_ts(1))

            @:ALLOCATE(pb_ts(2)%sf(idwbuff(1)%beg:idwbuff(1)%beg + 1, &
                idwbuff(2)%beg:idwbuff(2)%beg + 1, &
                idwbuff(3)%beg:idwbuff(3)%beg + 1, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(pb_ts(2))

            @:ALLOCATE(rhs_pb(idwbuff(1)%beg:idwbuff(1)%beg + 1, &
                idwbuff(2)%beg:idwbuff(2)%beg + 1, &
                idwbuff(3)%beg:idwbuff(3)%beg + 1, 1:nnode, 1:nb))
        else
            @:ALLOCATE(pb_ts(1)%sf(0,0,0,0,0))
            @:ACC_SETUP_SFs(pb_ts(1))

            @:ALLOCATE(pb_ts(2)%sf(0,0,0,0,0))
            @:ACC_SETUP_SFs(pb_ts(2))

            @:ALLOCATE(rhs_pb(0,0,0,0,0))
        end if

        @:ALLOCATE(mv_ts(1:2))

        if (qbmm .and. (.not. polytropic)) then
            @:ALLOCATE(mv_ts(1)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(mv_ts(1))

            @:ALLOCATE(mv_ts(2)%sf(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(mv_ts(2))

            @:ALLOCATE(rhs_mv(idwbuff(1)%beg:idwbuff(1)%end, &
                idwbuff(2)%beg:idwbuff(2)%end, &
                idwbuff(3)%beg:idwbuff(3)%end, 1:nnode, 1:nb))

        else if (qbmm .and. polytropic) then
            @:ALLOCATE(mv_ts(1)%sf(idwbuff(1)%beg:idwbuff(1)%beg + 1, &
                idwbuff(2)%beg:idwbuff(2)%beg + 1, &
                idwbuff(3)%beg:idwbuff(3)%beg + 1, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(mv_ts(1))

            @:ALLOCATE(mv_ts(2)%sf(idwbuff(1)%beg:idwbuff(1)%beg + 1, &
                idwbuff(2)%beg:idwbuff(2)%beg + 1, &
                idwbuff(3)%beg:idwbuff(3)%beg + 1, 1:nnode, 1:nb))
            @:ACC_SETUP_SFs(mv_ts(2))

            @:ALLOCATE(rhs_mv(idwbuff(1)%beg:idwbuff(1)%beg + 1, &
                idwbuff(2)%beg:idwbuff(2)%beg + 1, &
                idwbuff(3)%beg:idwbuff(3)%beg + 1, 1:nnode, 1:nb))
        else
            @:ALLOCATE(mv_ts(1)%sf(0,0,0,0,0))
            @:ACC_SETUP_SFs(mv_ts(1))

            @:ALLOCATE(mv_ts(2)%sf(0,0,0,0,0))
            @:ACC_SETUP_SFs(mv_ts(2))

            @:ALLOCATE(rhs_mv(0,0,0,0,0))
        end if

        ! Allocating the cell-average RHS variables
        @:ALLOCATE(rhs_vf(1:sys_size))
        @:PREFER_GPU(rhs_vf)

        if (igr) then
            do i = 1, sys_size
                @:ALLOCATE(rhs_vf(i)%sf(-1:m+1,-1:n+1,-1:p+1))
                @:ACC_SETUP_SFs(rhs_vf(i))
                @:PREFER_GPU(rhs_vf(i)%sf)
            end do
        else
            do i = 1, sys_size
                @:ALLOCATE(rhs_vf(i)%sf(0:m, 0:n, 0:p))
                @:ACC_SETUP_SFs(rhs_vf(i))
            end do
        end if

        ! Opening and writing the header of the run-time information file
        if (proc_rank == 0 .and. run_time_info) then
            call s_open_run_time_information_file()
        end if

        if (cfl_dt) then
            @:ALLOCATE(max_dt(0:m, 0:n, 0:p))
        end if

        ! Allocating arrays to store the bc types
        @:ALLOCATE(bc_type(1:num_dims,1:2))

        @:ALLOCATE(bc_type(1,1)%sf(0:0,0:n,0:p))
        @:ALLOCATE(bc_type(1,2)%sf(0:0,0:n,0:p))
        if (n > 0) then
            @:ALLOCATE(bc_type(2,1)%sf(-buff_size:m+buff_size,0:0,0:p))
            @:ALLOCATE(bc_type(2,2)%sf(-buff_size:m+buff_size,0:0,0:p))
            if (p > 0) then
                @:ALLOCATE(bc_type(3,1)%sf(-buff_size:m+buff_size,-buff_size:n+buff_size,0:0))
                @:ALLOCATE(bc_type(3,2)%sf(-buff_size:m+buff_size,-buff_size:n+buff_size,0:0))
            end if
        end if

        do i = 1, num_dims
            do j = 1, 2
                @:ACC_SETUP_SFs(bc_type(i,j))
            end do
        end do

        if (any(time_stepper == (/1, 2, 3, 4/))) then
            ! temporary array index for TVD RK
            if (time_stepper == 1) then
                stor = 1
            else
                stor = 2
            end if

            ! TVD RK coefficients
            @:ALLOCATE (rk_coef(time_stepper, 4))
            if (time_stepper == 1) then
                rk_coef(1, :) = (/1._wp, 0._wp, 1._wp, 1._wp/)
            else if (time_stepper == 2) then
                rk_coef(1, :) = (/1._wp, 0._wp, 1._wp, 1._wp/)
                rk_coef(2, :) = (/1._wp, 1._wp, 1._wp, 2._wp/)
            else if (time_stepper == 3) then
                rk_coef(1, :) = (/1._wp, 0._wp, 1._wp, 1._wp/)
                rk_coef(2, :) = (/1._wp, 3._wp, 1._wp, 4._wp/)
                rk_coef(3, :) = (/2._wp, 1._wp, 2._wp, 3._wp/)
            end if
            $:GPU_UPDATE(device='[rk_coef, stor]')
        end if

        @:ALLOCATE (rho_avg_x(-1:m, 0:n, 0:p), vel_avg_rms_x(-1:m, 0:n, 0:p), H_avg_x(-1:m, 0:n, 0:p), gamma_avg_x(-1:m, 0:n, 0:p))
        @:ALLOCATE (vel_avg_x(-1:m, 0:n, 0:p, 1:num_dims))
        @:ALLOCATE (rho_avg_y(0:m, -1:n, 0:p), vel_avg_rms_y(0:m, -1:n, 0:p), H_avg_y(0:m, -1:n, 0:p), gamma_avg_y(0:m, -1:n, 0:p))
        @:ALLOCATE (vel_avg_y(0:m, -1:n, 0:p, 1:num_dims))

    end subroutine s_initialize_time_steppers_module

    impure subroutine s_tvd_rk(t_step, time_avg, nstage)
#ifdef _CRAYFTN
        !DIR$ OPTIMIZE (-haggress)
#endif
        integer, intent(in) :: t_step
        real(wp), intent(inout) :: time_avg
        integer, intent(in) :: nstage

        integer :: i, j, k, l, q, s !< Generic loop iterator
        real(wp) :: start, finish
        integer :: dest

        call cpu_time(start)
        call nvtxStartRange("TIMESTEP")

        ! Adaptive dt: initial stage
        if (adap_dt) call s_adaptive_dt_bubble(1)

        do s = 1, nstage
            call s_compute_rhs(q_cons_ts(1)%vf, q_T_sf, q_prim_vf, bc_type, rhs_vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, t_step, time_avg, s)

            if (s == 1) then
                if (run_time_info) then
                    if (igr) then
                        call s_write_run_time_information(q_cons_ts(1)%vf, t_step)
                    else
                        call s_write_run_time_information(q_prim_vf, t_step)
                    end if
                end if

                if (probe_wrt) then
                    call s_time_step_cycling(t_step)
                end if

                if (cfl_dt) then
                    if (mytime >= t_stop) return
                else
                    if (t_step == t_step_stop) return
                end if
            end if

            if (bubbles_lagrange .and. .not. adap_dt) call s_update_lagrange_tdv_rk(stage=s)
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            if (s == 1 .and. nstage > 1) then
                                q_cons_ts(stor)%vf(i)%sf(j, k, l) = &
                                    q_cons_ts(1)%vf(i)%sf(j, k, l)
                            end if
                            if (igr) then
                                q_cons_ts(1)%vf(i)%sf(j, k, l) = &
                                    (rk_coef(s, 1)*q_cons_ts(1)%vf(i)%sf(j, k, l) &
                                     + rk_coef(s, 2)*q_cons_ts(stor)%vf(i)%sf(j, k, l) &
                                     + rk_coef(s, 3)*rhs_vf(i)%sf(j, k, l))/rk_coef(s, 4)
                            else
                                q_cons_ts(1)%vf(i)%sf(j, k, l) = &
                                    (rk_coef(s, 1)*q_cons_ts(1)%vf(i)%sf(j, k, l) &
                                     + rk_coef(s, 2)*q_cons_ts(stor)%vf(i)%sf(j, k, l) &
                                     + rk_coef(s, 3)*dt*rhs_vf(i)%sf(j, k, l))/rk_coef(s, 4)
                            end if
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            !Evolve pb and mv for non-polytropic qbmm
            if (qbmm .and. (.not. polytropic)) then
                $:GPU_PARALLEL_LOOP(collapse=5)
                do i = 1, nb
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                do q = 1, nnode
                                    if (s == 1 .and. nstage > 1) then
                                        pb_ts(stor)%sf(j, k, l, q, i) = &
                                            pb_ts(1)%sf(j, k, l, q, i)
                                        mv_ts(stor)%sf(j, k, l, q, i) = &
                                            mv_ts(1)%sf(j, k, l, q, i)
                                    end if
                                    pb_ts(1)%sf(j, k, l, q, i) = &
                                        (rk_coef(s, 1)*pb_ts(1)%sf(j, k, l, q, i) &
                                         + rk_coef(s, 2)*pb_ts(stor)%sf(j, k, l, q, i) &
                                         + rk_coef(s, 3)*dt*rhs_pb(j, k, l, q, i))/rk_coef(s, 4)
                                    mv_ts(1)%sf(j, k, l, q, i) = &
                                        (rk_coef(s, 1)*mv_ts(1)%sf(j, k, l, q, i) &
                                         + rk_coef(s, 2)*mv_ts(stor)%sf(j, k, l, q, i) &
                                         + rk_coef(s, 3)*dt*rhs_mv(j, k, l, q, i))/rk_coef(s, 4)
                                end do
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            if (bodyForces) call s_apply_bodyforces(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, rk_coef(s, 3)*dt/rk_coef(s, 4))

            if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(1)%vf)

            if (model_eqns == 3 .and. (.not. relax)) then
                call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)
            end if

            if (adv_n) call s_comp_alpha_from_n(q_cons_ts(1)%vf)

            if (ib) then
                ! check if any IBMS are moving, and if so, update the markers, ghost points, levelsets, and levelset norms
                if (moving_immersed_boundary_flag) then
                    do i = 1, num_ibs
                        if (s == 1) then
                            patch_ib(i)%step_vel = patch_ib(i)%vel
                            patch_ib(i)%step_angular_vel = patch_ib(i)%angular_vel
                            patch_ib(i)%step_angles = patch_ib(i)%angles
                            patch_ib(i)%step_x_centroid = patch_ib(i)%x_centroid
                            patch_ib(i)%step_y_centroid = patch_ib(i)%y_centroid
                            patch_ib(i)%step_z_centroid = patch_ib(i)%z_centroid
                        end if

                        if (patch_ib(i)%moving_ibm > 0) then
                            patch_ib(i)%vel = (rk_coef(s, 1)*patch_ib(i)%step_vel + rk_coef(s, 2)*patch_ib(i)%vel)/rk_coef(s, 4)
                            patch_ib(i)%angular_vel = (rk_coef(s, 1)*patch_ib(i)%step_angular_vel + rk_coef(s, 2)*patch_ib(i)%angular_vel)/rk_coef(s, 4)

                            if (patch_ib(i)%moving_ibm == 2) then ! if we are using two-way coupling, apply force and torque
                                ! compute the force and torque on the IB from the fluid
                                call s_compute_ib_forces(q_prim_vf(E_idx))

                                ! update the velocity from the force value
                                patch_ib(i)%vel = patch_ib(i)%vel + rk_coef(s, 3)*dt*(patch_ib(i)%force/patch_ib(i)%mass)/rk_coef(s, 4)

                                ! update the angular velocity with the torque value
                                patch_ib(i)%angular_vel = (patch_ib(i)%angular_vel*patch_ib(i)%moment) + (rk_coef(s, 3)*dt*patch_ib(i)%torque/rk_coef(s, 4)) ! add the torque to the angular momentum
                                call s_compute_moment_of_inertia(i, patch_ib(i)%angular_vel) ! update the moment of inertia to be based on the direction of the angular momentum
                                patch_ib(i)%angular_vel = patch_ib(i)%angular_vel/patch_ib(i)%moment ! convert back to angular velocity with the new moment of inertia
                            end if

                            ! Update the angle of the IB
                            patch_ib(i)%angles = (rk_coef(s, 1)*patch_ib(i)%step_angles + rk_coef(s, 2)*patch_ib(i)%angles + rk_coef(s, 3)*patch_ib(i)%angular_vel*dt)/rk_coef(s, 4)

                            ! Update the position of the IB
                            patch_ib(i)%x_centroid = (rk_coef(s, 1)*patch_ib(i)%step_x_centroid + rk_coef(s, 2)*patch_ib(i)%x_centroid + rk_coef(s, 3)*patch_ib(i)%vel(1)*dt)/rk_coef(s, 4)
                            patch_ib(i)%y_centroid = (rk_coef(s, 1)*patch_ib(i)%step_y_centroid + rk_coef(s, 2)*patch_ib(i)%y_centroid + rk_coef(s, 3)*patch_ib(i)%vel(2)*dt)/rk_coef(s, 4)
                            patch_ib(i)%z_centroid = (rk_coef(s, 1)*patch_ib(i)%step_z_centroid + rk_coef(s, 2)*patch_ib(i)%z_centroid + rk_coef(s, 3)*patch_ib(i)%vel(3)*dt)/rk_coef(s, 4)
                        end if
                    end do
                    call s_update_mib(num_ibs, levelset, levelset_norm)
                end if
                if (qbmm .and. .not. polytropic) then
                    call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf, pb_ts(1)%sf, mv_ts(1)%sf)
                else
                    call s_ibm_correct_state(q_cons_ts(1)%vf, q_prim_vf)
                end if
            end if
        end do

        ! Adaptive dt: final stage
        if (adap_dt) call s_adaptive_dt_bubble(3)

        call nvtxEndRange
        call cpu_time(finish)

        wall_time = abs(finish - start)

        if (t_step >= 2) then
            wall_time_avg = (wall_time + (t_step - 2)*wall_time_avg)/(t_step - 1)
        else
            wall_time_avg = 0._wp
        end if

    end subroutine s_tvd_rk

    subroutine s_dts(t_step, time_avg)
        integer, intent(in) :: t_step
        real(wp), intent(inout) :: time_avg
        real(wp) :: start, finish

        integer :: iter
        logical :: dts_conv
        real(wp) :: phi

        ! Start subroutine
        call cpu_time(start)
        call nvtxStartRange("TIMESTEP")

        ! Initialize q_cons_pts(1)
        call s_dts_initialize()

        ! Perform pseudo time iteration
        iter = 0; dts_conv = .false.
        do while (.true.)
            call s_dts_iteration()
            if (dts_conv) exit
            if (iter >= dts_iter_max) call s_mpi_abort("iter >= dts_iter_max")
        end do

        ! End subroutine
        call nvtxEndRange
        call cpu_time(finish)

        wall_time = abs(finish - start)

        if (t_step >= 2) then
            wall_time_avg = (wall_time + (t_step - 2)*wall_time_avg)/(t_step - 1)
        else
            wall_time_avg = 0._wp
        end if

    contains
        subroutine s_dts_initialize()
            integer :: i, j, k, l

            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = idwbuff(3)%beg, idwbuff(3)%end
                    do k = idwbuff(2)%beg, idwbuff(2)%end
                        do j = idwbuff(1)%beg, idwbuff(1)%end
                            q_cons_pts(1)%vf(i)%sf(j, k, l) = q_cons_ts(1)%vf(i)%sf(j, k, l)
                            q_cons_pts(2)%vf(i)%sf(j, k, l) = q_cons_ts(1)%vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do

            if (t_step == 0) then
                phi = 0._wp
            else
                phi = 0.5_wp
            end if
        end subroutine s_dts_initialize

        subroutine s_dts_iteration()
            ! Compute pseudo time variables
            call s_dts_update_pseudo_var_rk()

            ! Check convergence
            call s_dts_check_convergence()

            ! Update physical time variables
            if (dts_conv) call s_dts_update_physical_vars()
        end subroutine s_dts_iteration

        subroutine s_dts_update_pseudo_var_rk()
            real(wp), dimension(sys_size) :: dq
            real(wp), dimension(3) :: coeffs
            real(wp) :: dtp, dtp0
            integer :: nstage
            integer :: i, j, k, l, s

            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = idwbuff(3)%beg, idwbuff(3)%end
                    do k = idwbuff(2)%beg, idwbuff(2)%end
                        do j = idwbuff(1)%beg, idwbuff(1)%end
                            q_cons_pts(2)%vf(i)%sf(j, k, l) = q_cons_pts(1)%vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
            
            nstage = 3
            coeffs = (/0.1918_wp, 0.4929_wp, 1.0_wp/)

            do s = 1, nstage
                ! Compute RHS
                call s_compute_rhs(q_cons_pts(1)%vf, q_T_sf, q_prim_vf, bc_type, rhs_vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, t_step, time_avg, 1)

                $:GPU_PARALLEL_LOOP(collapse=3)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            ! Compute RHS for pseudo time step
                            do i = 1, sys_size
                                dq(i) = rhs_vf(i)%sf(j, k, l) &
                                        - ((1._wp + phi)*(q_cons_pts(1)%vf(i)%sf(j, k, l) - q_cons_ts(1)%vf(i)%sf(j, k, l)) &
                                           + phi*(q_cons_ts(1)%vf(i)%sf(j, k, l) - q_cons_ts(2)%vf(i)%sf(j, k, l)))/dt
                                ! NaN checker
                                if (dq(i) /= dq(i)) call s_mpi_abort("dq is NaN")
                            end do                                

                            ! Preconditioning (if flagged) and determining pseudo time step size
                            call s_dts_aux(dq, dtp, j, k, l)
                            if (s == 1) dtp0 = dtp

                            ! call s_dts_residual_smoothing()

                            ! Update pseudo time variables
                            do i = 1, sys_size
                                q_cons_pts(1)%vf(i)%sf(j, k, l) = q_cons_pts(2)%vf(i)%sf(j, k, l) + coeffs(s)*dtp0*dq(i)
                            end do
                        end do
                    end do
                end do
            end do
        end subroutine s_dts_update_pseudo_var_rk

        subroutine s_dts_update_pseudo_var_bdf2()
            real(wp), dimension(sys_size) :: dq
            real(wp) :: bb
            real(wp) :: phi, dtp
            integer :: i, j, k, l

            ! Compute RHS
            call s_compute_rhs(q_cons_pts(1)%vf, q_T_sf, q_prim_vf, bc_type, rhs_vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, t_step, time_avg, 1)

            if (t_step == 0) then
                phi = 0._wp
            else
                phi = 0.5_wp
            end if

            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        ! Compute RHS for pseudo time step
                        do i = 1, sys_size
                            dq(i) = rhs_vf(i)%sf(j, k, l) &
                                    - ((1._wp + phi)*(q_cons_pts(1)%vf(i)%sf(j, k, l) - q_cons_ts(1)%vf(i)%sf(j, k, l)) &
                                       + phi*(q_cons_ts(1)%vf(i)%sf(j, k, l) - q_cons_ts(2)%vf(i)%sf(j, k, l)))/dt

                            ! NaN checker
                            if (dq(i) /= dq(i)) then
                                print *, i, j, k, l, iter, rhs_vf(i)%sf(j, k, l), q_cons_pts(1)%vf(i)%sf(j, k, l), q_prim_vf(i)%sf(j, k, l)
                                call s_mpi_abort("dq is NaN")
                            end if
                        end do

                        ! Preconditioning (if flagged) and determining pseudo time step size
                        call s_dts_aux(dq, dtp, j, k, l)

                        ! Update pseudo time variables
                        do i = 1, sys_size
                            q_cons_pts(1)%vf(i)%sf(j, k, l) = q_cons_pts(1)%vf(i)%sf(j, k, l) + dtp*dq(i)
                        end do
                    end do
                end do
            end do
        end subroutine s_dts_update_pseudo_var_bdf2

        subroutine s_dts_check_convergence()
            real(wp) :: max_err, max_err_glb

            ! Compute error
            call s_dts_compute_error(max_err)

#ifdef MFC_MPI
            call s_mpi_allreduce_max(max_err, max_err_glb)
#else
            max_err_glb = max_err
#endif

            if (proc_rank == 0) write(99,*) iter, max_err_glb
            if (proc_rank == 0) print *, iter, max_err_glb

            ! Check convergence
            if (max_err_glb < 1._wp) then
                if (proc_rank == 0) print *, "converged at pseudo-time iteration: ", iter
                dts_conv = .true.
            else
                iter = iter + 1
            end if
        end subroutine s_dts_check_convergence

        subroutine s_dts_compute_error(max_err)
            real(wp), intent(out) :: max_err
            real(wp) :: max_err_tmp
            real(wp), dimension(sys_size) :: err
            integer :: i, j, k, l

            call s_compute_rhs(q_cons_pts(1)%vf, q_T_sf, q_prim_vf, bc_type, rhs_vf, pb_ts(1)%sf, rhs_pb, mv_ts(1)%sf, rhs_mv, t_step, time_avg, 1)

            max_err = 0._wp
            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        ! Compute error
                        max_err_tmp = 0._wp
                        do i = 1, sys_size
                            err(i) = rhs_vf(i)%sf(j, k, l)*dt &
                                     - ((1._wp + phi)*(q_cons_pts(1)%vf(i)%sf(j, k, l) - q_cons_ts(1)%vf(i)%sf(j, k, l)) &
                                        + phi*(q_cons_ts(1)%vf(i)%sf(j, k, l) - q_cons_ts(2)%vf(i)%sf(j, k, l)))
                            ! max_err_tmp = max_err_tmp + abs(err(i))/(dts_a_tol)
                            max_err_tmp = max_err_tmp + abs(err(i))/(dts_a_tol + dts_r_tol*abs(q_cons_ts(1)%vf(i)%sf(j, k, l)))
                        end do
                        ! Update max err
                        if (max_err_tmp > max_err) max_err = max_err_tmp
                    end do
                end do
            end do
        end subroutine s_dts_compute_error

        subroutine s_dts_update_physical_vars()
            integer :: i, j, k, l
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_cons_ts(2)%vf(i)%sf(j, k, l) = q_cons_ts(1)%vf(i)%sf(j, k, l)
                            q_cons_ts(1)%vf(i)%sf(j, k, l) = q_cons_pts(1)%vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
        end subroutine s_dts_update_physical_vars

        subroutine s_dts_aux(dq, dtp, j, k, l)
            real(wp), dimension(sys_size), intent(inout) :: dq
            real(wp), intent(inout) :: dtp
            real(wp), dimension(sys_size) :: pcond
            real(wp) :: rho        !< Cell-avg. density
            real(wp), dimension(num_vels) :: vel        !< Cell-avg. velocity
            real(wp) :: vel_sum    !< Cell-avg. velocity sum
            real(wp) :: pres       !< Cell-avg. pressure
            real(wp), dimension(num_fluids) :: alpha      !< Cell-avg. volume fraction
            real(wp) :: gamma      !< Cell-avg. sp. heat ratio
            real(wp) :: pi_inf     !< Cell-avg. liquid stiffness function
            real(wp) :: qv         !< Cell-avg. fluid reference energy
            real(wp) :: c          !< Cell-avg. sound speed
            real(wp) :: H          !< Cell-avg. enthalpy
            real(wp), dimension(2) :: Re         !< Cell-avg. Reynolds numbers
            real(wp) :: rho_K, gamma_K, pi_inf_K, vel_sum_K
            real(wp) :: beta, lambda
            real(wp) :: ds
            real(wp) :: aa, bb, cc
            integer :: i, j, k, l

            call s_compute_enthalpy(q_prim_vf, pres, rho, gamma, pi_inf, Re, H, alpha, vel, vel_sum, qv, j, k, l)
            call s_compute_speed_of_sound(pres, rho, gamma, pi_inf, H, alpha, vel_sum, 0._wp, c, qv)
            beta = min(max(dts_cutoff**2._wp, vel_sum/c**2._wp), 1._wp)

            ! Largest eigenvalue
            if (preconditioning) then
                lambda = 0.5_wp*((1._wp + beta)*sqrt(vel_sum) + sqrt((1._wp - beta)**2._wp*vel_sum + 4._wp*beta*c**2._wp))
            else
                lambda = sqrt(vel_sum) + c
            end if

            ! Pseudo time step size
            ds = minval(dx); if (n > 0) ds = min(ds, minval(dy)); if (p > 0) ds = min(ds, minval(dz))
            dtp = dts_cfl/(lambda/ds)

            aa = 1.5_wp*dtp/dt; bb = 1._wp/(1._wp + aa); cc = bb/(1._wp + aa*beta)
            if (preconditioning) then
                ! Continuity
                do i = contxb, contxe
                    pcond(i) = (beta - 1._wp)*vel_sum/2._wp*cc
                end do

                ! Momentum
                do i = 1, num_dims
                    pcond(momxb - 1 + i) = (1._wp - beta)*vel(i)*cc
                end do

                ! Energy
                pcond(E_idx) = beta*cc/bb

                ! Volume fractions
                do i = 1, num_fluids
                    pcond(advxb - 1 + i) = (1._wp - beta)*(gammas(i)*pres + pi_infs(i))*cc
                end do
            else
                pcond = 1._wp
                pcond(E_idx) = bb
            end if

            do i = 1, sys_size
                if (i == E_idx) then
                    dq(E_idx) = sum(pcond*dq)
                else
                    dq(i) = dq(i)*bb
                end if
            end do
        end subroutine s_dts_aux

        subroutine s_dts_residual_smoothing()
            ! integer :: id

            ! do id = 1, num_dims
            !     ! Compute Jacobian
            !     call s_dts_compute_flux_jacobian(&
            !       !----->In
            !       q_cons_left, &
            !       q_cons_right, &
            !       !----->Out
            !       fjacobian_p, &
            !       fjacobian_m, &
            !       id &
            !     )
            ! end do

            ! ! Compute smoothed residual
            ! call s_dts_red_black_solver(&
            !   !----->In
            !   fjacobian_p_total, &
            !   fjacobian_m_total, &
            !   !----->Out
            !   q_smoothed_residual &
            ! )
        end subroutine s_dts_residual_smoothing

        subroutine s_dts_compute_flux_jacobian()
            real(wp), dimension(sys_size, sys_size) :: fjacobian

            ! ! Compute Jacobian - pseudo code now
            ! fjacobian = 0._wp
            ! fjacobian(1, 2) = 1._wp
            ! fjacobian(2, 1) = -vel_avg(1)**2._wp / rho_avg + vel_avg_rms/(2._wp*rho_avg*gamma_avg)
            ! fjacobian(2, 2) = (1._wp - 1._wp / gamma_avg) * vel_avg(1) / rho_avg
            ! fjacobian(2, 3) = -vel_avg(2) / (rho_avg * gamma_avg)
            ! fjacobian(2, 4) = 
            ! fjacobian(3, 1) = -vel_avg(1)*vel_avg(2) / rho_avg
            ! fjacobian(3, 3) = vel_avg(1) / rho_avg
            ! fjacobian(4, 1) = 
            ! fjacobian(4, 2) = 
            ! fjacobian(4, 3) = 
            ! fjacobian(4, 4) = 

            ! fjacobian_p = 0.5_wp*(fjacobian + abs(fjacobian))
            ! fjacobian_m = 0.5_wp*(fjacobian - abs(fjacobian))

            ! ! Compute LHS contribution
            ! call s_dts_add_flux_jacobian(&
            !   !----->In
            !   fjacobian_p, &
            !   !----->Inout
            !   fjacobian_p_total &
            ! )

            ! ! Compute LHS contribution
            ! call s_dts_add_flux_jacobian(&
            !   !----->In
            !   fjacobian_m, &
            !   !----->Inout
            !   fjacobian_m_total &
            ! )            
        end subroutine s_dts_compute_flux_jacobian

        subroutine s_dts_add_flux_jacobian()

        end subroutine s_dts_add_flux_jacobian

        subroutine s_dts_red_black_solver()

        end subroutine s_dts_red_black_solver

    end subroutine s_dts

    ! subroutine s_dts_update_states(rho_avg, vel_avg, vel_avg_rms, H_avg, gamma_avg, 
    !                                                           j, k, l, norm_dir)
    !     if (norm_dir == 1) then
    !         rho_avg_x(j, k, l) = rho_avg
    !         do i = 1, num_dims
    !             vel_avg_x(j, k, l, i) = vel_avg(i)
    !         end do
    !         vel_avg_rms_x(j, k, l) = vel_avg_rms
    !         H_avg_x(j, k, l) = H_avg
    !         gamma_avg_x(j, k, l) = gamma_avg
    !       else if (norm_dir == 2) then
    !         rho_avg_y(j, k, l) = rho_avg
    !         do i = 1, num_dims
    !             vel_avg_y(j, k, l, i) = vel_avg(i)
    !         end do
    !         vel_avg_rms_y(j, k, l) = vel_avg_rms
    !         H_avg_y(j, k, l) = H_avg
    !         gamma_avg_y(j, k, l) = gamma_avg
    !       end if
    ! end subroutine s_dts_update_states

    !> Bubble source part in Strang operator splitting scheme
        !! @param t_step Current time-step
    impure subroutine s_adaptive_dt_bubble(stage)

        integer, intent(in) :: stage

        type(vector_field) :: gm_alpha_qp

        call s_convert_conservative_to_primitive_variables( &
            q_cons_ts(1)%vf, &
            q_T_sf, &
            q_prim_vf, &
            idwint)

        if (bubbles_euler) then

            call s_compute_bubble_EE_source(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, divu)
            call s_comp_alpha_from_n(q_cons_ts(1)%vf)

        elseif (bubbles_lagrange) then

            call s_populate_variables_buffers(bc_type, q_prim_vf, pb_ts(1)%sf, mv_ts(1)%sf)
            call s_compute_bubble_EL_dynamics(q_prim_vf, stage)
            call s_transfer_data_to_tmp()
            call s_smear_voidfraction()
            if (stage == 3) then
                if (lag_params%write_bubbles_stats) call s_calculate_lag_bubble_stats()
                if (lag_params%write_bubbles) then
                    $:GPU_UPDATE(host='[gas_p,gas_mv,intfc_rad,intfc_vel]')
                    call s_write_lag_particles(mytime)
                end if
                call s_write_void_evol(mytime)
            end if

        end if

    end subroutine s_adaptive_dt_bubble

    impure subroutine s_compute_dt()

        real(wp) :: rho        !< Cell-avg. density
        real(wp), dimension(num_vels) :: vel        !< Cell-avg. velocity
        real(wp) :: vel_sum    !< Cell-avg. velocity sum
        real(wp) :: pres       !< Cell-avg. pressure
        real(wp), dimension(num_fluids) :: alpha      !< Cell-avg. volume fraction
        real(wp) :: gamma      !< Cell-avg. sp. heat ratio
        real(wp) :: pi_inf     !< Cell-avg. liquid stiffness function
        real(wp) :: qv         !< Cell-avg. fluid reference energy
        real(wp) :: c          !< Cell-avg. sound speed
        real(wp) :: H          !< Cell-avg. enthalpy
        real(wp), dimension(2) :: Re         !< Cell-avg. Reynolds numbers
        type(vector_field) :: gm_alpha_qp

        real(wp) :: dt_local
        integer :: j, k, l !< Generic loop iterators

        if (.not. igr) then
            call s_convert_conservative_to_primitive_variables( &
                q_cons_ts(1)%vf, &
                q_T_sf, &
                q_prim_vf, &
                idwint)
        end if

        $:GPU_PARALLEL_LOOP(collapse=3, private='[vel, alpha, Re, rho, vel_sum, pres, gamma, pi_inf, c, H, qv]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    if (igr) then
                        call s_compute_enthalpy(q_cons_ts(1)%vf, pres, rho, gamma, pi_inf, Re, H, alpha, vel, vel_sum, qv, j, k, l)
                    else
                        call s_compute_enthalpy(q_prim_vf, pres, rho, gamma, pi_inf, Re, H, alpha, vel, vel_sum, qv, j, k, l)
                    end if

                    ! Compute mixture sound speed
                    call s_compute_speed_of_sound(pres, rho, gamma, pi_inf, H, alpha, vel_sum, 0._wp, c, qv)

                    call s_compute_dt_from_cfl(vel, c, max_dt, rho, Re, j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        #:call GPU_PARALLEL(copyout='[dt_local]', copyin='[max_dt]')
            dt_local = minval(max_dt)
        #:endcall GPU_PARALLEL

        if (num_procs == 1) then
            dt = dt_local
        else
            call s_mpi_allreduce_min(dt_local, dt)
        end if

        $:GPU_UPDATE(device='[dt]')

    end subroutine s_compute_dt

    !> This subroutine applies the body forces source term at each
        !! Runge-Kutta stage
    subroutine s_apply_bodyforces(q_cons_vf, q_prim_vf_in, rhs_vf_in, ldt)

        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(1:sys_size), intent(in) :: q_prim_vf_in
        type(scalar_field), dimension(1:sys_size), intent(inout) :: rhs_vf_in

        real(wp), intent(in) :: ldt !< local dt

        integer :: i, j, k, l

        call nvtxStartRange("RHS-BODYFORCES")
        call s_compute_body_forces_rhs(q_prim_vf_in, q_cons_vf, rhs_vf_in)

        $:GPU_PARALLEL_LOOP(collapse=4)
        do i = momxb, E_idx
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l) + &
                                                   ldt*rhs_vf_in(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        call nvtxEndRange

    end subroutine s_apply_bodyforces

    !> This subroutine saves the temporary q_prim_vf vector
        !!      into the q_prim_ts vector that is then used in p_main
        !! @param t_step current time-step
    subroutine s_time_step_cycling(t_step)

        integer, intent(in) :: t_step

        integer :: i, j, k, l !< Generic loop iterator

        if (t_step == t_step_start) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_prim_ts2(2)%vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        elseif (t_step == t_step_start + 1) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_prim_ts2(1)%vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        elseif (t_step == t_step_start + 2) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_prim_ts1(2)%vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        elseif (t_step == t_step_start + 3) then
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_prim_ts1(1)%vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        else ! All other timesteps
            $:GPU_PARALLEL_LOOP(collapse=4)
            do i = 1, sys_size
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_prim_ts2(2)%vf(i)%sf(j, k, l) = q_prim_ts2(1)%vf(i)%sf(j, k, l)
                            q_prim_ts2(1)%vf(i)%sf(j, k, l) = q_prim_ts1(2)%vf(i)%sf(j, k, l)
                            q_prim_ts1(2)%vf(i)%sf(j, k, l) = q_prim_ts1(1)%vf(i)%sf(j, k, l)
                            q_prim_ts1(1)%vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_time_step_cycling

    !> Module deallocation and/or disassociation procedures
    impure subroutine s_finalize_time_steppers_module
#ifdef FRONTIER_UNIFIED
        use hipfort
        use hipfort_hipmalloc
        use hipfort_check
#endif
        integer :: i, j !< Generic loop iterators

        ! Deallocating the cell-average conservative variables
#if defined(__NVCOMPILER_GPU_UNIFIED_MEM)
        do j = 1, sys_size
            @:DEALLOCATE(q_cons_ts(1)%vf(j)%sf)
            @:DEALLOCATE(q_cons_pts(1)%vf(j)%sf)
            if (num_ts == 2) then
                if (nv_uvm_out_of_core) then
                    nullify (q_cons_ts(2)%vf(j)%sf)
                    nullify (q_cons_pts(2)%vf(j)%sf)
                else
                    @:DEALLOCATE(q_cons_ts(2)%vf(j)%sf)
                    @:DEALLOCATE(q_cons_pts(2)%vf(j)%sf)
                end if
            end if
        end do
        if (num_ts == 2 .and. nv_uvm_out_of_core) then
            deallocate (q_cons_ts_pool_host)
            deallocate (q_cons_pts_pool_host)
        end if
#elif defined(FRONTIER_UNIFIED)
        do i = 1, num_ts
            do j = 1, sys_size
                nullify (q_cons_ts(i)%vf(j)%sf)
                nullify (q_cons_pts(i)%vf(j)%sf)
            end do
        end do
#ifdef MFC_MIXED_PRECISION
        call hipCheck(hipHostFree_(c_loc(q_cons_ts_pool_host)))
        nullify (q_cons_ts_pool_host)
        call hipCheck(hipFree_(c_loc(q_cons_ts_pool_device)))
        nullify (q_cons_ts_pool_device)
#else
        call hipCheck(hipHostFree(q_cons_ts_pool_host))
        call hipCheck(hipFree(q_cons_ts_pool_device))
        call hipCheck(hipHostFree(q_cons_pts_pool_host))
        call hipCheck(hipFree(q_cons_pts_pool_device))
#endif
#else
        do i = 1, num_ts
            do j = 1, sys_size
                @:DEALLOCATE(q_cons_ts(i)%vf(j)%sf)
                @:DEALLOCATE(q_cons_pts(i)%vf(j)%sf)
            end do
        end do
#endif
        do i = 1, num_ts
            @:DEALLOCATE(q_cons_ts(i)%vf)
            @:DEALLOCATE(q_cons_pts(i)%vf)
        end do

        @:DEALLOCATE(q_cons_ts)
        @:DEALLOCATE(q_cons_pts)

        ! Deallocating the cell-average primitive ts variables
        if (probe_wrt) then
            do i = 1, num_probe_ts
                do j = 1, sys_size
                    @:DEALLOCATE(q_prim_ts1(i)%vf(j)%sf,q_prim_ts2(i)%vf(j)%sf )
                end do
                @:DEALLOCATE(q_prim_ts1(i)%vf, q_prim_ts2(i)%vf)
            end do
            @:DEALLOCATE(q_prim_ts1, q_prim_ts2)
        end if

        if (.not. igr) then
            ! Deallocating the cell-average primitive variables
            do i = 1, adv_idx%end
                @:DEALLOCATE(q_prim_vf(i)%sf)
            end do

            if (mhd) then
                do i = B_idx%beg, B_idx%end
                    @:DEALLOCATE(q_prim_vf(i)%sf)
                end do
            end if

            if (elasticity) then
                do i = stress_idx%beg, stress_idx%end
                    @:DEALLOCATE(q_prim_vf(i)%sf)
                end do
            end if

            if (hyperelasticity) then
                do i = xibeg, xiend + 1
                    @:DEALLOCATE(q_prim_vf(i)%sf)
                end do
            end if

            if (cont_damage) then
                @:DEALLOCATE(q_prim_vf(damage_idx)%sf)
            end if

            if (bubbles_euler) then
                do i = bub_idx%beg, bub_idx%end
                    @:DEALLOCATE(q_prim_vf(i)%sf)
                end do
            end if

            if (model_eqns == 3) then
                do i = internalEnergies_idx%beg, internalEnergies_idx%end
                    @:DEALLOCATE(q_prim_vf(i)%sf)
                end do
            end if
        end if

        @:DEALLOCATE(q_prim_vf)

        ! Deallocating the cell-average RHS variables
        do i = 1, sys_size
            @:DEALLOCATE(rhs_vf(i)%sf)
        end do

        @:DEALLOCATE(rhs_vf)

        ! Writing the footer of and closing the run-time information file
        if (proc_rank == 0 .and. run_time_info) then
            call s_close_run_time_information_file()
        end if

    end subroutine s_finalize_time_steppers_module

end module m_time_steppers
