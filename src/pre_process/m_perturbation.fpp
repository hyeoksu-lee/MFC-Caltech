!>
!! @file m_perturbation.fpp
!! @brief Contains module m_perturbation

!> @brief This module contains subroutines that compute perturbations to the
!!              initial mean flow fields.
module m_perturbation

    use m_derived_types         ! Definitions of the derived types

    use m_global_parameters     ! Global parameters for the code

    use m_mpi_proxy              !< Message passing interface (MPI) module proxy

    use m_eigen_solver          ! Subroutines to solve eigenvalue problem for
    ! complex general matrix

    use ieee_arithmetic

    implicit none

    integer :: mixlayer_nvar ! Number of variables in linear stability analysis solver for mixing layer
    integer, allocatable, dimension(:) :: mixlayer_var ! Index of variables in linear stability analysis solver
    integer :: nbpm, nbp ! Number of grid cell boundary points in y-direction
    integer :: n0, nbpm0, nbp0 ! Number of grid cell boundary points in y-direction
    real(wp) :: dy0
    real(wp), allocatable, dimension(:) :: y_cb0, dy_cb0
    integer :: mixlayer_bc_fd ! Order of finite difference applied at the boundaries of mixing layer
    integer :: n_bc_skip ! Number of points skipped in the linear stability analysis due to the boundary condition

contains

    subroutine s_initialize_perturbation_module()

        if (mixlayer_perturb) then
            mixlayer_bc_fd = 2
            n0 = 255
            nbpm0 = n0 + 1
            nbp0 = n0 + 2
            nbpm = n + 1
            nbp = n + 2
            allocate (y_cb0(0:nbpm0))
            allocate (dy_cb0(0:nbpm0 - 1))
            
            if (model_eqns == 2 .and. num_fluids == 1) then
                n_bc_skip = mixlayer_bc_fd*2
                mixlayer_nvar = 5 ! 1 continuity + 3 momentum + 1 energy
                allocate (mixlayer_var(mixlayer_nvar))

                mixlayer_var(1) = contxb
                mixlayer_var(2) = momxb
                mixlayer_var(3) = momxb + 1
                mixlayer_var(4) = momxb + 2
                mixlayer_var(5) = momxb + 3
            end if
        end if

    end subroutine s_initialize_perturbation_module

    subroutine s_perturb_sphere(q_prim_vf)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        integer :: i, j, k, l !< generic loop operators

        real(wp) :: perturb_alpha
        real(wp) :: alpha_unadv
        real(wp) :: rand_real
        call random_seed()

        do k = 0, p
            do j = 0, n
                do i = 0, m
                    call random_number(rand_real)

                    perturb_alpha = q_prim_vf(E_idx + perturb_sph_fluid)%sf(i, j, k)

                    ! Perturb partial density fields to match perturbed volume fraction fields
                    !    IF ((perturb_alpha >= 25e-2_wp) .AND. (perturb_alpha <= 75e-2_wp)) THEN
                    if ((perturb_alpha /= 0._wp) .and. (perturb_alpha /= 1._wp)) then

                        ! Derive new partial densities
                        do l = 1, num_fluids
                            q_prim_vf(l)%sf(i, j, k) = q_prim_vf(E_idx + l)%sf(i, j, k)*fluid_rho(l)
                        end do

                    end if
                end do
            end do
        end do

    end subroutine s_perturb_sphere

    subroutine s_perturb_surrounding_flow(q_prim_vf)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        integer :: i, j, k, l !<  generic loop iterators

        real(wp) :: perturb_alpha
        real(wp) :: rand_real
        call random_seed()

        ! Perturb partial density or velocity of surrounding flow by some random small amount of noise
        do k = 0, p
            do j = 0, n
                do i = 0, m
                    perturb_alpha = q_prim_vf(E_idx + perturb_flow_fluid)%sf(i, j, k)
                    call random_number(rand_real)
                    rand_real = rand_real*perturb_flow_mag
                    q_prim_vf(mom_idx%beg)%sf(i, j, k) = (1._wp + rand_real)*q_prim_vf(mom_idx%beg)%sf(i, j, k)
                    q_prim_vf(mom_idx%end)%sf(i, j, k) = rand_real*q_prim_vf(mom_idx%beg)%sf(i, j, k)
                    if (bubbles_euler) then
                        q_prim_vf(alf_idx)%sf(i, j, k) = (1._wp + rand_real)*q_prim_vf(alf_idx)%sf(i, j, k)
                    end if
                end do
            end do
        end do

    end subroutine s_perturb_surrounding_flow

    subroutine s_perturb_random_noise(q_prim_vf)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        integer :: i, j, k, l !<  generic loop iterators
        real(wp) :: uratio
        real(wp) :: rand_real

        uratio = 1._wp/patch_icpp(1)%vel(1)

        call random_seed()

        do k = 0, p
            do j = 0, n
                do i = 0, m
                    call random_number(rand_real)
                    q_prim_vf(momxb)%sf(i, j, k) = q_prim_vf(momxb)%sf(i, j, k) + (2._wp * rand_real - 1._wp) * 0.01_wp / uratio! u
                    call random_number(rand_real)
                    q_prim_vf(momxb + 1)%sf(i, j, k) = q_prim_vf(momxb + 1)%sf(i, j, k) + (2._wp * rand_real - 1._wp) * 0.01_wp / uratio ! v
                    if (p > 0) then
                        call random_number(rand_real)
                        q_prim_vf(momxb + 2)%sf(i, j, k) = q_prim_vf(momxb + 2)%sf(i, j, k) + (2._wp * rand_real - 1._wp) * 0.01_wp / uratio ! w
                    end if
                    call random_number(rand_real)
                    q_prim_vf(E_idx)%sf(i, j, k) = q_prim_vf(E_idx)%sf(i, j, k) + (2._wp * rand_real - 1._wp) * 0.01_wp / uratio**2._wp  ! p
                end do
            end do
        end do

    end subroutine s_perturb_random_noise

    !>  This subroutine computes velocity perturbations for a temporal mixing
        !!              layer with hypertangent mean streamwise velocity profile
        !!              obtained from linear stability analysis. For a 2D case,
        !!              instability waves with spatial wavenumbers, (4,0), (2,0),
        !!              and (1,0) are superposed. For a 3D waves, (4,4), (4,-4),
        !!              (2,2), (2,-2), (1,1), (1,-1) areadded on top of 2D waves.
    subroutine s_superposition_instability_wave(q_prim_vf)
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(wp), dimension(mixlayer_nvar, 0:m, 0:n, 0:p) :: wave, wave1, wave2, wave_tmp
        real(wp) :: uratio, Ldomain
        real(wp) :: nR3bar
        real(wp), dimension(6) :: shift
        integer :: i, j, k, q

        uratio = 1._wp/patch_icpp(1)%vel(1)
        Ldomain = mixlayer_domain
        
        ! Generate base grid
        call s_generate_base_grid()

        ! Generate waves
        wave = 0._wp
        wave1 = 0._wp
        wave2 = 0._wp

        ! Compute 2D waves
        call s_instability_wave(2*pi*4.0_wp/Ldomain, 0._wp, wave_tmp, 0._wp)
        print *, "1"
        wave1 = wave1 + wave_tmp
        call s_instability_wave(2*pi*2.0_wp/Ldomain, 0._wp, wave_tmp, 0._wp)
        print *, "2"
        wave1 = wave1 + wave_tmp
        call s_instability_wave(2*pi*1.0_wp/Ldomain, 0._wp, wave_tmp, 0._wp)
        print *, "3"
        wave1 = wave1 + wave_tmp
        wave = wave1*0.05_wp

        if (mixlayer_shift == 1) then
            shift(1) = 2*pi*11._wp/31._wp; shift(2) = 2*pi*13._wp/31._wp; shift(3) = 2*pi*17._wp/31._wp;
            shift(4) = 2*pi*19._wp/31._wp; shift(5) = 2*pi*23._wp/31._wp; shift(6) = 2*pi*29._wp/31._wp;
            print *, "shift 1", (shift(i), i=1,6)
        else if (mixlayer_shift == 2) then
            shift(1) = 2*pi*7._wp/61._wp;  shift(2) = 2*pi*11._wp/61._wp; shift(3) = 2*pi*19._wp/61._wp;
            shift(4) = 2*pi*41._wp/61._wp; shift(5) = 2*pi*53._wp/61._wp; shift(6) = 2*pi*59._wp/61._wp;
            print *, "shift 2", (shift(i), i=1,6)
        else if (mixlayer_shift == 3) then
            shift(1) = 2*pi*17._wp/53._wp; shift(2) = 2*pi*19._wp/53._wp; shift(3) = 2*pi*31._wp/53._wp;
            shift(4) = 2*pi*47._wp/53._wp; shift(5) = 2*pi*29._wp/53._wp; shift(6) = 2*pi*3._wp/53._wp;
            print *, "shift 3", (shift(i), i=1,6)
        else if (mixlayer_shift == 4) then
            shift(1) = 2*pi*13._wp/43._wp; shift(2) = 2*pi*11._wp/43._wp; shift(3) = 2*pi*39._wp/43._wp;
            shift(4) = 2*pi*29._wp/43._wp; shift(5) = 2*pi*23._wp/43._wp; shift(6) = 2*pi*19._wp/43._wp;
            print *, "shift 4", (shift(i), i=1,6)
        else if (mixlayer_shift == 5) then
            shift(1) = 2*pi*19._wp/37._wp; shift(2) = 2*pi*31._wp/37._wp; shift(3) = 2*pi*29._wp/37._wp;
            shift(4) = 2*pi*3._wp/37._wp;  shift(5) = 2*pi*23._wp/37._wp; shift(6) = 2*pi*11._wp/37._wp;
            print *, "shift 5", (shift(i), i=1,6)
        else if (mixlayer_shift == 6) then
            shift(1) = 2*pi*2._wp/67._wp;  shift(2) = 2*pi*53._wp/67._wp; shift(3) = 2*pi*29._wp/67._wp;
            shift(4) = 2*pi*31._wp/67._wp; shift(5) = 2*pi*13._wp/67._wp; shift(6) = 2*pi*17._wp/67._wp;
            print *, "shift 6", (shift(i), i=1,6)
        else if (mixlayer_shift == 7) then
            shift(1) = 2*pi*0.1377641919_wp; shift(2) = 2*pi*0.6152866056_wp; shift(3) = 2*pi*0.0511069403_wp;
            shift(4) = 2*pi*0.3927737151_wp; shift(5) = 2*pi*0.0006762054_wp; shift(6) = 2*pi*0.1675698069_wp;
            print *, "shift 7", (shift(i), i=1,6)
        else if (mixlayer_shift == 8) then
            shift(1) = 2*pi*0.8556119067_wp; shift(2) = 2*pi*0.9513815031_wp; shift(3) = 2*pi*0.4328811743_wp;
            shift(4) = 2*pi*0.4420256142_wp; shift(5) = 2*pi*0.7374684184_wp; shift(6) = 2*pi*0.2477236891_wp;
            print *, "shift 8", (shift(i), i=1,6)
        else if (mixlayer_shift == 9) then
            shift(1) = 2*pi*0.5576752220_wp; shift(2) = 2*pi*0.8501618422_wp; shift(3) = 2*pi*0.3463564289_wp;
            shift(4) = 2*pi*0.4973941163_wp; shift(5) = 2*pi*0.9409184387_wp; shift(6) = 2*pi*0.9994186173_wp;
            print *, "shift 9", (shift(i), i=1,6)
        else if (mixlayer_shift == 10) then
            shift(1) = 2*pi*0.4899695672_wp; shift(2) = 2*pi*0.8505344420_wp; shift(3) = 2*pi*0.3721308791_wp;
            shift(4) = 2*pi*0.7329680488_wp; shift(5) = 2*pi*0.2234294319_wp; shift(6) = 2*pi*0.9827279516_wp;
            print *, "shift 10", (shift(i), i=1,6)
        else if (mixlayer_shift == 11) then
            shift(1) = 2*pi*0.73936086_wp; shift(2) = 2*pi*0.53321271_wp; shift(3) = 2*pi*0.07294918_wp;
            shift(4) = 2*pi*0.54818918_wp; shift(5) = 2*pi*0.86098990_wp; shift(6) = 2*pi*0.54304788_wp;
            print *, "shift 11", (shift(i), i=1,6)
        else if (mixlayer_shift == 12) then
            shift(1) = 2*pi*0.63102723_wp; shift(2) = 2*pi*0.55273117_wp; shift(3) = 2*pi*0.67750635_wp;
            shift(4) = 2*pi*0.88108276_wp; shift(5) = 2*pi*0.06290314_wp; shift(6) = 2*pi*0.60971874_wp;
            print *, "shift 12", (shift(i), i=1,6)
        else if (mixlayer_shift == 13) then
            shift(1) = 2*pi*0.28108231_wp; shift(2) = 2*pi*0.57440735_wp; shift(3) = 2*pi*0.11556867_wp;
            shift(4) = 2*pi*0.20123427_wp; shift(5) = 2*pi*0.54014151_wp; shift(6) = 2*pi*0.68002622_wp;
            print *, "shift 13", (shift(i), i=1,6)
        else if (mixlayer_shift == 14) then
            shift(1) = 2*pi*0.97913344_wp; shift(2) = 2*pi*0.73720290_wp; shift(3) = 2*pi*0.96491401_wp;
            shift(4) = 2*pi*0.19457755_wp; shift(5) = 2*pi*0.38104128_wp; shift(6) = 2*pi*0.42362568_wp;
            print *, "shift 14", (shift(i), i=1,6)
        else if (mixlayer_shift == 15) then
            shift(1) = 2*pi*0.82113501_wp; shift(2) = 2*pi*0.70251862_wp; shift(3) = 2*pi*0.34801304_wp;
            shift(4) = 2*pi*0.47004825_wp; shift(5) = 2*pi*0.94827641_wp; shift(6) = 2*pi*0.03544434_wp;
            print *, "shift 15", (shift(i), i=1,6)
        else if (mixlayer_shift == 16) then
            shift(1) = 2*pi*0.04775845_wp; shift(2) = 2*pi*0.21268208_wp; shift(3) = 2*pi*0.78296638_wp;
            shift(4) = 2*pi*0.42997254_wp; shift(5) = 2*pi*0.47059390_wp; shift(6) = 2*pi*0.80020549_wp;
            print *, "shift 16", (shift(i), i=1,6)
        else if (mixlayer_shift == 17) then
            shift(1) = 2*pi*0.21239576_wp; shift(2) = 2*pi*0.85857176_wp; shift(3) = 2*pi*0.46040151_wp;
            shift(4) = 2*pi*0.64376148_wp; shift(5) = 2*pi*0.95176286_wp; shift(6) = 2*pi*0.36655479_wp;
            print *, "shift 17", (shift(i), i=1,6)
        else if (mixlayer_shift == 18) then
            shift(1) = 2*pi*0.85103539_wp; shift(2) = 2*pi*0.64516243_wp; shift(3) = 2*pi*0.52306323_wp;
            shift(4) = 2*pi*0.98260012_wp; shift(5) = 2*pi*0.03886104_wp; shift(6) = 2*pi*0.11332287_wp;
            print *, "shift 18", (shift(i), i=1,6)
        else if (mixlayer_shift == 19) then
            shift(1) = 2*pi*0.79283785_wp; shift(2) = 2*pi*0.47085308_wp; shift(3) = 2*pi*0.19101602_wp;
            shift(4) = 2*pi*0.96079179_wp; shift(5) = 2*pi*0.56967998_wp; shift(6) = 2*pi*0.73977715_wp;
            print *, "shift 19", (shift(i), i=1,6)
        else if (mixlayer_shift == 20) then
            shift(1) = 2*pi*0.62122631_wp; shift(2) = 2*pi*0.50381444_wp; shift(3) = 2*pi*0.34548855_wp;
            shift(4) = 2*pi*0.24801029_wp; shift(5) = 2*pi*0.08947161_wp; shift(6) = 2*pi*0.27224491_wp;
            print *, "shift 20", (shift(i), i=1,6)
        end if

        if (p > 0) then
            ! Compute 3D waves with phase shifts.
            call s_instability_wave(2*pi*4.0_wp/Ldomain, 2*pi*4.0_wp/Ldomain, wave_tmp, shift(1))
            print *, "4"
            wave2 = wave2 + wave_tmp
            call s_instability_wave(2*pi*2.0_wp/Ldomain, 2*pi*2.0_wp/Ldomain, wave_tmp, shift(2))
            print *, "5"
            wave2 = wave2 + wave_tmp
            call s_instability_wave(2*pi*1.0_wp/Ldomain, 2*pi*1.0_wp/Ldomain, wave_tmp, shift(3))
            print *, "6"
            wave2 = wave2 + wave_tmp
            call s_instability_wave(2*pi*4.0_wp/Ldomain, -2*pi*4.0_wp/Ldomain, wave_tmp, shift(4))
            print *, "7"
            wave2 = wave2 + wave_tmp
            call s_instability_wave(2*pi*2.0_wp/Ldomain, -2*pi*2.0_wp/Ldomain, wave_tmp, shift(5))
            print *, "8"
            wave2 = wave2 + wave_tmp
            call s_instability_wave(2*pi*1.0_wp/Ldomain, -2*pi*1.0_wp/Ldomain, wave_tmp, shift(6))
            print *, "9"
            wave2 = wave2 + wave_tmp
            wave = wave + 0.15_wp*wave2
        end if

        ! Superpose velocity perturbuations (instability waves) to the velocity field
        do k = 0, p
            do j = 0, n
                do i = 0, m
                    q_prim_vf(contxb)%sf(i, j, k) = q_prim_vf(contxb)%sf(i, j, k) + wave(mixlayer_var(1), i, j, k) ! rho
                    q_prim_vf(momxb)%sf(i, j, k) = q_prim_vf(momxb)%sf(i, j, k) + wave(mixlayer_var(2), i, j, k)/uratio ! u
                    q_prim_vf(momxb + 1)%sf(i, j, k) = q_prim_vf(momxb + 1)%sf(i, j, k) + wave(mixlayer_var(3), i, j, k)/uratio ! v
                    if (p > 0) then
                        q_prim_vf(momxb + 2)%sf(i, j, k) = q_prim_vf(momxb + 2)%sf(i, j, k) + wave(mixlayer_var(4), i, j, k)/uratio ! w
                    end if
                    q_prim_vf(E_idx)%sf(i, j, k) = q_prim_vf(E_idx)%sf(i, j, k) + wave(mixlayer_var(5), i, j, k)/uratio**2._wp ! p

                    if (bubbles_euler .and. (.not. qbmm)) then
                        nR3bar = 0._wp
                        do q = 1, nb
                            call s_compute_equilibrium_state(q_prim_vf(E_idx)%sf(i, j, k), R0(q), q_prim_vf(bub_idx%rs(q))%sf(i, j, k))
                            nR3bar = nR3bar + weight(q)*(q_prim_vf(n_idx)%sf(i, j, k)*q_prim_vf(bub_idx%rs(q))%sf(i, j, k))**3._wp
                        end do
                        if (.not. decouple) then
                            q_prim_vf(alf_idx)%sf(i, j, k) = (4._wp*pi*nR3bar)/(3._wp*q_prim_vf(n_idx)%sf(i, j, k)**2._wp)
                        end if
                    end if
                end do
            end do
        end do

    end subroutine s_superposition_instability_wave

    !>  This subroutine computes equilibrium bubble radius of the perturbed pressure field
    subroutine s_compute_equilibrium_state(fP, fR0, fR)
        real(wp), intent(in) :: fP, fR0
        real(wp), intent(inout) :: fR
        real(wp) :: f0, f1
        real(wp) :: gam_b
        integer :: ii, jj

        gam_b = 1._wp + 1._wp/fluid_pp(num_fluids + 1)%gamma

        ! Loop
        ii = 1
        do while (.true.)

            f0 = (Ca + 2._wp/(Web*fR0))*(fR0/fR)**(3._wp*gam_b) - 2._wp/(Web*fR) + 1._wp - Ca - fP
            f1 = -3._wp*gam_b*(Ca + 2._wp/(Web*fR0))*(fR0/fR)**(3._wp*gam_b + 1._wp) / fR0 + 2._wp/(Web*fR**2._wp)

            if (abs(f0) <= 1e-10_wp) then
                ! Converged
                exit
            else
                ! Update radius
                fR = fR - f0/f1
            end if

            ! Failed case
            if (ieee_is_nan(f0) .or. &
                ieee_is_nan(f1) .or. &
                ii > 1000 .or. &
                fR < 0._wp) then

                print *, "Failed to compute equilibrium radius"
                print *, ii, f0, f1, fR0, fR, gam_b, Ca, Web, fP
                call s_mpi_abort()
                fR = fR0
                exit
            end if

            ii = ii + 1
        end do

    end subroutine s_compute_equilibrium_state

    !>  This subroutine computes instability waves for a given set of spatial
        !!              wavenumbers (alpha, beta) in x and z directions.
        !!              The eigenvalue problem is derived from the linearized
        !!              Euler equations with parallel mean flow assumption
        !!              (See Sandham 1989 PhD thesis for details).
    subroutine s_instability_wave(alpha, beta, wave, shift)
        real(wp), intent(in) :: alpha, beta !<  spatial wavenumbers
        real(wp), dimension(mixlayer_nvar, 0:m, 0:n, 0:p), intent(inout) :: wave !< instability wave
        real(wp), intent(in) :: shift !< phase shift
        real(wp), dimension(0:nbpm0) :: u_mean !<  mean density and velocity profiles
        real(wp) :: rho_mean, p_mean !< mean density and pressure
        real(wp), dimension(0:nbpm0, 0:nbpm0) :: d !< differential operator in y dir
        real(wp) :: gam, pi_inf, mach, c1, adv
        real(wp) :: xratio, uratio
        integer :: i, j !<  generic loop iterators

        xratio = mixlayer_vel_coef
        uratio = 1._wp/patch_icpp(1)%vel(1)

        ! Set fluid flow properties
        if (bubbles_euler) then
            adv = patch_icpp(1)%alpha(num_fluids)
        else
            adv = 0._wp
        end if
        gam = 1._wp + 1._wp/fluid_pp(1)%gamma
        pi_inf = fluid_pp(1)%pi_inf*(gam - 1._wp)/gam*uratio**2
        rho_mean = patch_icpp(1)%alpha_rho(1)
        p_mean = patch_icpp(1)%pres*uratio**2
        c1 = sqrt((gam*(p_mean + pi_inf))/(rho_mean*(1._wp - adv)))
        mach = 1._wp/c1

        ! Assign mean profiles
        do j = 0, nbpm0
            u_mean(j) = tanh(y_cb0(j)*xratio)
        end do

        ! Compute differential operator in y-dir
        ! based on 2nd order central difference
        d = 0._wp
        d(0, 0) = -1._wp/((y_cb0(1) - y_cb0(0))*xratio)
        d(0, 1) = 1._wp/((y_cb0(1) - y_cb0(0))*xratio)
        do j = 1, nbpm0 - 1
            d(j, j - 1) = -1._wp/((y_cb0(j + 1) - y_cb0(j - 1))*xratio)
            d(j, j + 1) = 1._wp/((y_cb0(j + 1) - y_cb0(j - 1))*xratio)
        end do
        d(nbpm0, nbpm0 - 1) = -1._wp/((y_cb0(nbpm0) - y_cb0(nbpm0 - 1))*xratio)
        d(nbpm0, nbpm0) = 1._wp/((y_cb0(nbpm0) - y_cb0(nbpm0 - 1))*xratio)

        ! Compute
        call s_solve_linear_system(alpha, beta, u_mean, rho_mean, p_mean, d, gam, pi_inf, mach, wave, shift)

    end subroutine s_instability_wave

    !> This subroutine solves linear system from linear stability analysis and
        !!              generate instability waves for the given set of spatial
        !!              wave numbers and phase shift.
    subroutine s_solve_linear_system(alpha, beta, u_mean, rho_mean, p_mean, d, gam, pi_inf, mach, wave, shift)
        real(wp), intent(in) :: alpha, beta !<  spatial wavenumbers
        real(wp), dimension(0:nbpm0), intent(in) :: u_mean !<  mean velocity profiles
        real(wp), intent(in) :: rho_mean, p_mean !< mean density and pressure
        real(wp), dimension(0:nbpm0, 0:nbpm0), intent(in) :: d !< differential operator in y dir
        real(wp), intent(in) :: gam, pi_inf, mach, shift
        real(wp), dimension(mixlayer_nvar, 0:m, 0:n, 0:p), intent(inout) :: wave

        real(wp), dimension(0:nbpm0) :: drho_mean, du_mean !< y-derivatives of mean profiles
        real(wp), dimension(0:mixlayer_nvar*nbp0 - 1, 0:mixlayer_nvar*nbp0 - 1) :: ar, ai    !< matrices for eigenvalue problem
        real(wp), dimension(0:mixlayer_nvar*nbp0 - 1, 0:mixlayer_nvar*nbp0 - 1) :: br, bi, ci !< matrices for eigenvalue problem
        real(wp), dimension(0:mixlayer_nvar*n0 - n_bc_skip - 1, 0:mixlayer_nvar*n0 - n_bc_skip - 1) :: hr, hi    !< matrices for eigenvalue problem

        real(wp), dimension(0:mixlayer_nvar*n0 - n_bc_skip - 1, 0:mixlayer_nvar*n0 - n_bc_skip - 1) :: zr, zi !< eigenvectors
        real(wp), dimension(0:mixlayer_nvar*n0 - n_bc_skip - 1) :: wr, wi !< eigenvalues
        real(wp), dimension(0:mixlayer_nvar*n0 - n_bc_skip - 1) :: fv1, fv2, fv3 !< temporary memory

        integer :: ierr
        integer :: i, j, k, l !<  generic loop iterators
        integer :: ii, jj !< block matrix indices

        ! Compute y-derivatives of rho and u
        do j = 0, nbpm0
            drho_mean(j) = 0
            du_mean(j) = 0
            do k = 0, nbpm0
                drho_mean(j) = 0._wp
                du_mean(j) = du_mean(j) + d(j, k)*u_mean(k)
            end do
        end do

        ! Compute B and C, then A = B + C. Here, A is the matrix for the linear
        ! systems of equation (i.e. we are going to solve x for Ax = lambda x).
        ! Here, B includes components of A without differential operator, and
        ! C includes components of A with differential operator.
        br = 0._wp
        bi = 0._wp
        ci = 0._wp
        do j = 0, nbpm0
            ii = mixlayer_var(1); jj = mixlayer_var(1); br((ii - 1)*nbp0 + j, (jj - 1)*nbp0 + j) = alpha*u_mean(j); 
            ii = mixlayer_var(1); jj = mixlayer_var(2); br((ii - 1)*nbp0 + j, (jj - 1)*nbp0 + j) = alpha*rho_mean; 
            ii = mixlayer_var(1); jj = mixlayer_var(3); bi((ii - 1)*nbp0 + j, (jj - 1)*nbp0 + j) = -drho_mean(j); 
            ii = mixlayer_var(1); jj = mixlayer_var(4); br((ii - 1)*nbp0 + j, (jj - 1)*nbp0 + j) = beta*rho_mean; 
            ii = mixlayer_var(2); jj = mixlayer_var(2); br((ii - 1)*nbp0 + j, (jj - 1)*nbp0 + j) = alpha*u_mean(j); 
            ii = mixlayer_var(2); jj = mixlayer_var(3); bi((ii - 1)*nbp0 + j, (jj - 1)*nbp0 + j) = -du_mean(j); 
            ii = mixlayer_var(2); jj = mixlayer_var(5); br((ii - 1)*nbp0 + j, (jj - 1)*nbp0 + j) = alpha/rho_mean; 
            ii = mixlayer_var(3); jj = mixlayer_var(3); br((ii - 1)*nbp0 + j, (jj - 1)*nbp0 + j) = alpha*u_mean(j); 
            ii = mixlayer_var(4); jj = mixlayer_var(4); br((ii - 1)*nbp0 + j, (jj - 1)*nbp0 + j) = alpha*u_mean(j); 
            ii = mixlayer_var(4); jj = mixlayer_var(5); br((ii - 1)*nbp0 + j, (jj - 1)*nbp0 + j) = beta/rho_mean; 
            ii = mixlayer_var(5); jj = mixlayer_var(2); br((ii - 1)*nbp0 + j, (jj - 1)*nbp0 + j) = gam*(p_mean + pi_inf)*alpha; 
            ii = mixlayer_var(5); jj = mixlayer_var(4); br((ii - 1)*nbp0 + j, (jj - 1)*nbp0 + j) = gam*(p_mean + pi_inf)*beta; 
            ii = mixlayer_var(5); jj = mixlayer_var(5); br((ii - 1)*nbp0 + j, (jj - 1)*nbp0 + j) = alpha*u_mean(j); 
            do k = 0, nbpm0
                ii = mixlayer_var(1); jj = mixlayer_var(3); ci((ii - 1)*nbp0 + j, (jj - 1)*nbp0 + k) = -rho_mean*d(j, k); 
                ii = mixlayer_var(3); jj = mixlayer_var(5); ci((ii - 1)*nbp0 + j, (jj - 1)*nbp0 + k) = -d(j, k)/rho_mean; 
                ii = mixlayer_var(5); jj = mixlayer_var(3); ci((ii - 1)*nbp0 + j, (jj - 1)*nbp0 + k) = -gam*(p_mean + pi_inf)*d(j, k); 
            end do
        end do
        ar = br
        ai = bi + ci

        ! Apply BC to ar and ai matrices
        if (bc_y%beg == -6 .and. bc_y%end == -6) then
            ! Nonreflecting subsonic buffer BC
            call s_instability_nonreflecting_subsonic_buffer_bc(ar, ai, hr, hi, rho_mean, mach)
        end if

        ! Compute eigenvalues and eigenvectors
        call cg(mixlayer_nvar*n0 - n_bc_skip, mixlayer_nvar*n0 - n_bc_skip, hr, hi, wr, wi, zr, zi, fv1, fv2, fv3, ierr)

        ! Generate instability wave
        call s_generate_wave(wr, wi, zr, zi, rho_mean, mach, alpha, beta, wave, shift)

    end subroutine s_solve_linear_system

    !> This subroutine applies non-reflecting subsonic buffer boundary condition
        !!              to the linear system of equations (i.e. matrix A).
    subroutine s_instability_nonreflecting_subsonic_buffer_bc(ar, ai, hr, hi, rho_mean, mach)
        real(wp), dimension(0:mixlayer_nvar*nbp0 - 1, 0:mixlayer_nvar*nbp0 - 1), intent(inout) :: ar, ai    !< matrices for eigenvalue problem
        real(wp), dimension(0:mixlayer_nvar*n0 - n_bc_skip - 1, 0:mixlayer_nvar*n0 - n_bc_skip - 1), intent(out) :: hr, hi    !< matrices for eigenvalue problem
        real(wp), intent(in) :: rho_mean !<  mean density profiles
        real(wp), intent(in) :: mach
        real(wp), dimension(0:mixlayer_nvar*n0 - 1, 0:mixlayer_nvar*n0 - 1) :: fr, fi    !< matrices for eigenvalue problem
        real(wp), dimension(0:mixlayer_nvar*n0 - n_bc_skip - 1, 0:mixlayer_nvar*n0 - 1) :: gr, gi    !< matrices for eigenvalue problem
        integer :: i, j, k, l, ii, jj

        ! Condition 1: v = 0 at BC - no action required here

        ! Condition 2: du/dy = 0 at BC
        do j = 0, mixlayer_nvar*nbp0 - 1
            ! at y_domain%beg
            ii = mixlayer_var(1)*nbp0
            ar(j, ii + 1) = ar(j, ii + 1) + ar(j, ii)
            ai(j, ii + 1) = ai(j, ii + 1) + ai(j, ii)
            ! at y_domain%end
            ii = mixlayer_var(1)*nbp0 + nbp0 - 1
            ar(j, ii - 1) = ar(j, ii - 1) + ar(j, ii)
            ai(j, ii - 1) = ai(j, ii - 1) + ai(j, ii)
        end do

        ! Condition 3: dw/dy = 0 at BC
        do j = 0, mixlayer_nvar*nbp0 - 1
            ! at y_domain%beg
            ii = (mixlayer_var(3))*nbp0
            ar(j, ii + 1) = ar(j, ii + 1) + ar(j, ii)
            ai(j, ii + 1) = ai(j, ii + 1) + ai(j, ii)
            ! at y_domain%end
            ii = (mixlayer_var(3))*nbp0 + nbp0 - 1
            ar(j, ii - 1) = ar(j, ii - 1) + ar(j, ii)
            ai(j, ii - 1) = ai(j, ii - 1) + ai(j, ii)
        end do

        ! Condition 4: dp/dy +- rho c dv/dy = 0 at BC
        do j = 0, mixlayer_nvar*nbp0 - 1
            ! at y_domain%beg
            ii = mixlayer_var(4)*nbp0
            ar(j, ii + 1) = ar(j, ii + 1) + ar(j, ii)
            ai(j, ii + 1) = ai(j, ii + 1) + ai(j, ii)
            jj = mixlayer_var(2)*nbp0
            ar(j, jj + 1) = ar(j, jj + 1) + ar(j, ii)*rho_mean/mach
            ai(j, jj + 1) = ai(j, jj + 1) + ai(j, ii)*rho_mean/mach
            ! at y_domain%end
            ii = mixlayer_var(4)*nbp0 + nbp0 - 1
            ar(j, ii - 1) = ar(j, ii - 1) + ar(j, ii)
            ai(j, ii - 1) = ai(j, ii - 1) + ai(j, ii)
            jj = mixlayer_var(2)*nbp0 + nbp0 - 1
            ar(j, jj - 1) = ar(j, jj - 1) - ar(j, ii)*rho_mean/mach
            ai(j, jj - 1) = ai(j, jj - 1) - ai(j, ii)*rho_mean/mach
        end do

        ! Condition 5: c^2 drho/dy +- dp/dy = 0 at BC
        do j = 0, mixlayer_nvar*nbp0 - 1
            ! at y_domain%beg
            ii = 0
            ar(j, ii + 1) = ar(j, ii + 1) + ar(j, ii)
            ai(j, ii + 1) = ai(j, ii + 1) + ai(j, ii)
            jj = mixlayer_var(2)*nbp0
            ar(j, jj + 1) = ar(j, jj + 1) + ar(j, ii)*rho_mean*mach
            ai(j, jj + 1) = ai(j, jj + 1) + ai(j, ii)*rho_mean*mach
            ! at y_domain%end
            ii = nbp0 - 1
            ar(j, ii - 1) = ar(j, ii - 1) + ar(j, ii)
            ai(j, ii - 1) = ai(j, ii - 1) + ai(j, ii)
            jj = mixlayer_var(2)*nbp0 + nbp0 - 1
            ar(j, jj - 1) = ar(j, jj - 1) - ar(j, ii)*rho_mean*mach
            ai(j, jj - 1) = ai(j, jj - 1) - ai(j, ii)*rho_mean*mach
        end do

        ! Remove unnecessary rows of the matrix A (rho, u, v, w, p at the boundaries)
        fr = 0._wp
        fi = 0._wp
        do ii = 1, mixlayer_nvar
            do jj = 1, mixlayer_nvar
                do k = 0, n0 - 1
                    do l = 0, n0 - 1
                        fr((ii - 1)*n0 + k, (jj - 1)*n0 + l) = ar((ii - 1)*nbp0 + k + 1, (jj - 1)*nbp0 + l + 1)
                        fi((ii - 1)*n0 + k, (jj - 1)*n0 + l) = ai((ii - 1)*nbp0 + k + 1, (jj - 1)*nbp0 + l + 1)
                    end do
                end do
            end do
        end do

        gr = 0._wp
        gi = 0._wp
        do ii = 1, mixlayer_nvar
            do j = 0, mixlayer_nvar*n0 - 1
                if (ii <= mixlayer_var(2)) then
                    do k = 0, n0 - 1
                        gr((ii - 1)*n0 + k, j) = fr((ii - 1)*n0 + k, j)
                        gi((ii - 1)*n0 + k, j) = fi((ii - 1)*n0 + k, j)
                    end do
                elseif (ii == mixlayer_var(3)) then
                    do k = 0, n0 - n_bc_skip - 1
                        gr((ii - 1)*n0 + k, j) = fr((ii - 1)*n0 + k + mixlayer_bc_fd, j)
                        gi((ii - 1)*n0 + k, j) = fi((ii - 1)*n0 + k + mixlayer_bc_fd, j)
                    end do
                else
                    do k = 0, n0 - 1
                        gr((ii - 1)*n0 - n_bc_skip + k, j) = fr((ii - 1)*n0 + k, j)
                        gi((ii - 1)*n0 - n_bc_skip + k, j) = fi((ii - 1)*n0 + k, j)
                    end do
                end if
            end do
        end do

        hr = 0._wp
        hi = 0._wp
        do i = 0, mixlayer_nvar*n0 - n_bc_skip - 1
            do jj = 1, mixlayer_nvar
                if (jj <= mixlayer_var(2)) then
                    do k = 0, n0 - 1
                        hr(i, (jj - 1)*n0 + k) = gr(i, (jj - 1)*n0 + k)
                        hi(i, (jj - 1)*n0 + k) = gi(i, (jj - 1)*n0 + k)
                    end do
                elseif (jj == mixlayer_var(3)) then
                    do k = 0, n0 - n_bc_skip - 1
                        hr(i, (jj - 1)*n0 + k) = gr(i, (jj - 1)*n0 + k + mixlayer_bc_fd)
                        hi(i, (jj - 1)*n0 + k) = gi(i, (jj - 1)*n0 + k + mixlayer_bc_fd)
                    end do
                else
                    do k = 0, n0 - 1
                        hr(i, (jj - 1)*n0 - n_bc_skip + k) = gr(i, (jj - 1)*n0 + k)
                        hi(i, (jj - 1)*n0 - n_bc_skip + k) = gi(i, (jj - 1)*n0 + k)
                    end do
                end if
            end do
        end do

    end subroutine s_instability_nonreflecting_subsonic_buffer_bc

    !>  This subroutine generates an instability wave using the most unstable
        !!              eigenvalue and corresponding eigenvector among the
        !!              given set of eigenvalues and eigenvectors.
    subroutine s_generate_wave(wr, wi, zr, zi, rho_mean, mach, alpha, beta, wave, shift)
        real(wp), dimension(0:mixlayer_nvar*n0 - n_bc_skip - 1), intent(in) :: wr, wi !< eigenvalues
        real(wp), dimension(0:mixlayer_nvar*n0 - n_bc_skip - 1, 0:mixlayer_nvar*n0 - n_bc_skip - 1), intent(in) :: zr, zi !< eigenvectors
        real(wp), intent(in) :: rho_mean
        real(wp), dimension(mixlayer_nvar, 0:m, 0:n, 0:p), intent(inout) :: wave
        real(wp), intent(in) :: alpha, beta, mach, shift
        real(wp), dimension(0:mixlayer_nvar*n0 - n_bc_skip - 1) :: vr, vi, vnr, vni !< most unstable eigenvector
        real(wp), dimension(0:mixlayer_nvar*nbp0 - 1) :: xbr, xbi !< eigenvectors
        real(wp), dimension(0:mixlayer_nvar*nbp - 1) :: ubr, ubi
        real(wp), dimension(0:mixlayer_nvar*nbpm - 1) :: ucr, uci
        real(wp) :: ang, norm
        real(wp) :: tr, ti, cr, ci !< temporary memory
        real(wp) :: xratio
        integer :: idx
        integer :: i, j, k

        xratio = mixlayer_vel_coef

        ! Find the most unstable eigenvalue and corresponding eigenvector
        k = 0
        do i = 1, mixlayer_nvar*n0 - n_bc_skip - 1
            if (wi(i) > wi(k)) then
                k = i
            end if
        end do
        vr = zr(:, k)
        vi = zi(:, k)

        ! Normalize the eigenvector by its component with the largest modulus.
        norm = 0._wp
        do i = 0, mixlayer_nvar*n0 - n_bc_skip - 1
            if (sqrt(vr(i)**2 + vi(i)**2) > norm) then
                idx = i
                norm = sqrt(vr(i)**2 + vi(i)**2)
            end if
        end do

        tr = vr(idx)
        ti = vi(idx)
        do i = 0, mixlayer_nvar*n0 - n_bc_skip - 1
            call cdiv(vr(i), vi(i), tr, ti, cr, ci)
            vnr(i) = cr
            vni(i) = ci
        end do

        ! Reassign missing values at boundaries based on the boundary condition
        xbr = 0._wp
        xbi = 0._wp
        do i = 1, mixlayer_nvar
            if (i <= mixlayer_var(2)) then
                do k = 0, n0 - 1
                    xbr((i - 1)*nbp0 + k + 1) = vnr((i - 1)*n0 + k)
                    xbi((i - 1)*nbp0 + k + 1) = vni((i - 1)*n0 + k)
                end do
            elseif (i == mixlayer_var(3)) then
                do k = 0, n0 - n_bc_skip - 1
                    xbr((i - 1)*nbp0 + mixlayer_bc_fd + k + 1) = vnr((i - 1)*n0 + k)
                    xbi((i - 1)*nbp0 + mixlayer_bc_fd + k + 1) = vni((i - 1)*n0 + k)
                end do
            else
                do k = 0, n0 - 1
                    xbr((i - 1)*nbp0 + k + 1) = vnr((i - 1)*n0 - n_bc_skip + k)
                    xbi((i - 1)*nbp0 + k + 1) = vni((i - 1)*n0 - n_bc_skip + k)
                end do
            end if
        end do

        ! rho at boundaries
        xbr(0) = xbr(1) + xbr(mixlayer_var(2)*nbp0 + 1)*rho_mean*mach
        xbi(0) = xbi(1) + xbi(mixlayer_var(2)*nbp0 + 1)*rho_mean*mach
        xbr(nbp0 - 1) = xbr(n0) - xbr(mixlayer_var(2)*nbp0 + n0)*rho_mean*mach
        xbi(nbp0 - 1) = xbi(n0) - xbi(mixlayer_var(2)*nbp0 + n0)*rho_mean*mach

        ! u at boundaries
        xbr(mixlayer_var(1)*nbp0) = xbr(mixlayer_var(1)*nbp0 + 1)
        xbi(mixlayer_var(1)*nbp0) = xbi(mixlayer_var(1)*nbp0 + 1)
        xbr(mixlayer_var(1)*nbp0 + nbp0 - 1) = xbr(mixlayer_var(1)*nbp0 + n0)
        xbi(mixlayer_var(1)*nbp0 + nbp0 - 1) = xbi(mixlayer_var(1)*nbp0 + n0)

        ! w at boundaries
        xbr((mixlayer_var(3))*nbp0 + 0) = xbr((mixlayer_var(3))*nbp0 + 1)
        xbi((mixlayer_var(3))*nbp0 + 0) = xbi((mixlayer_var(3))*nbp0 + 1)
        xbr((mixlayer_var(3))*nbp0 + nbp0 - 1) = xbr((mixlayer_var(3))*nbp0 + n0)
        xbi((mixlayer_var(3))*nbp0 + nbp0 - 1) = xbi((mixlayer_var(3))*nbp0 + n0)

        ! p at boundaries
        xbr(mixlayer_var(4)*nbp0 + 0) = xbr(mixlayer_var(4)*nbp0 + 1) + xbr(mixlayer_var(2)*nbp0 + 1)*rho_mean/mach
        xbi(mixlayer_var(4)*nbp0 + 0) = xbi(mixlayer_var(4)*nbp0 + 1) + xbi(mixlayer_var(2)*nbp0 + 1)*rho_mean/mach
        xbr(mixlayer_var(4)*nbp0 + nbp0 - 1) = xbr(mixlayer_var(4)*nbp0 + n0) - xbr(mixlayer_var(2)*nbp0 + n0)*rho_mean/mach
        xbi(mixlayer_var(4)*nbp0 + nbp0 - 1) = xbi(mixlayer_var(4)*nbp0 + n0) - xbi(mixlayer_var(2)*nbp0 + n0)*rho_mean/mach

        ! Interpolate xbr and xbi to get ubr and ubi
        ! At boundaries
        do k = 0, 4
            ubr(k*nbp + 0) = xbr(k*nbp0 + 0)
            ubi(k*nbp + 0) = xbi(k*nbp0 + 0)
            ubr(k*nbp + nbpm) = xbr(k*nbp0 + nbpm0)
            ubi(k*nbp + nbpm) = xbi(k*nbp0 + nbpm0)
        end do

        do i = 1, n
            do j = 1, nbp0 - 1
                if (y_cb(i - 1) .ge. y_cb0(j - 1) .and. y_cb(i - 1) .lt. y_cb0(j) ) then
                    do k = 0, 4
                        ubr(k*nbp + i) = xbr(k*nbp0 + j - 1) + (xbr(k*nbp0 + j) - xbr(k*nbp0 + j - 1))/(y_cb0(j) - y_cb0(j - 1))*(y_cb(i - 1) - y_cb0(j - 1))
                        ubi(k*nbp + i) = xbi(k*nbp0 + j - 1) + (xbi(k*nbp0 + j) - xbi(k*nbp0 + j - 1))/(y_cb0(j) - y_cb0(j - 1))*(y_cb(i - 1) - y_cb0(j - 1))
                    end do
                    exit
                end if
            end do
        end do

        do i = 0, nbpm
            write(97,*) alpha, beta, i, y_cb(i - 1), ubr(i)
        end do

        do i = 0, nbpm0
            write(98,*) alpha, beta, i, y_cb0(i), ubr(i)
        end do

        ! Compute average to get cell-centered values
        ucr = 0._wp
        uci = 0._wp
        do i = 1, mixlayer_nvar
            do k = 0, n
                ucr((i - 1)*nbpm + k) = 5e-1_wp*(ubr((i - 1)*nbp + k) + ubr((i - 1)*nbp + k + 1))
                uci((i - 1)*nbpm + k) = 5e-1_wp*(ubi((i - 1)*nbp + k) + ubi((i - 1)*nbp + k + 1))
            end do
        end do

        ! Generate instability waves in x- and z-directions with phase shifts
        ! wave = Re(eigfunc * exp(i*(alpha*x + beta*z)))
        do i = 0, m
            do j = 0, n
                do k = 0, p
                    if (beta == 0) then
                        ang = alpha*(x_cc(i)*xratio)
                    else
                        ang = alpha*(x_cc(i)*xratio) + beta*(z_cc(k)*xratio) + shift
                    end if
                    wave(mixlayer_var(1), i, j, k) = ucr(j)*cos(ang) - uci(j)*sin(ang) ! rho
                    wave(mixlayer_var(2), i, j, k) = ucr(mixlayer_var(1)*nbpm + j)*cos(ang) - uci(mixlayer_var(1)*nbpm + j)*sin(ang) ! u
                    wave(mixlayer_var(3), i, j, k) = ucr(mixlayer_var(2)*nbpm + j)*cos(ang) - uci(mixlayer_var(2)*nbpm + j)*sin(ang) ! v
                    wave(mixlayer_var(4), i, j, k) = ucr(mixlayer_var(3)*nbpm + j)*cos(ang) - uci(mixlayer_var(3)*nbpm + j)*sin(ang) ! w
                    wave(mixlayer_var(5), i, j, k) = ucr(mixlayer_var(4)*nbpm + j)*cos(ang) - uci(mixlayer_var(4)*nbpm + j)*sin(ang) ! p
                end do
            end do
        end do

    end subroutine s_generate_wave

    subroutine s_generate_base_grid()

        real(wp) :: length  !< domain lengths
        real(wp) :: dy_tmp, y_a0, y_b0, a_y0
        integer :: i, j     !< generic loop operators

        length = abs(y_cb(n) - y_cb(-1))
        dy_tmp = length/real(nbpm0, wp)

        print *, y_cb(n), y_cb(-1), length

        do i = 0, nbpm0
            y_cb0(i) = y_cb(-1) + dy_tmp*real(i, wp)
        end do

        y_cb0 = y_cb0/length
        y_a0 = -5._wp/length
        y_b0 =  5._wp/length
        a_y0 = 5._wp

        do j = 1, 3
            do i = 0, nbpm0
                y_cb0(i) = y_cb0(i)/a_y0* &
                            (a_y0 + log(cosh(a_y0*(y_cb0(i) - y_a0))) &
                            + log(cosh(a_y0*(y_cb0(i) - y_b0))) &
                            - 2._wp*log(cosh(a_y0*(y_b0 - y_a0)/2._wp)))
            end do
        end do

        y_cb0 = y_cb0/(y_cb0(nbpm0) - y_cb0(0))*length
        dy_cb0 = y_cb0(1:nbpm0) - y_cb0(0:nbpm0 - 1)

        print *, y_cb0(nbpm0), y_cb0(0), y_cb0(nbpm0) - y_cb0(0)

    end subroutine s_generate_base_grid

end module m_perturbation
