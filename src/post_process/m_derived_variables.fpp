!>
!! @file m_derived_variables.f90
!! @brief Contains module m_derived_variables

!> @brief This module features subroutines that allow for the derivation of
!!      numerous flow variables from the conservative and primitive ones.
!!      Currently, the available derived variables include the unadvected
!!      volume fraction, specific heat ratio, liquid stiffness, speed of
!!      sound, vorticity and the numerical Schlieren function.

module m_derived_variables

    use m_derived_types         !< Definitions of the derived types

    use m_global_parameters     !< Global parameters for the code

    use m_mpi_proxy             !< Message passing interface (MPI) module proxy

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_variables_conversion

    implicit none

    private; public :: s_initialize_derived_variables_module, &
 s_derive_specific_heat_ratio, &
 s_derive_liquid_stiffness, &
 s_derive_sound_speed, &
 s_derive_flux_limiter, &
 s_derive_vorticity_component, &
 s_derive_qm, &
 s_derive_liutex, &
 s_apply_gaussian_filter, &
 s_add_paddings_xy, &
 s_add_paddings_real, &
 s_add_paddings_logical, &
 s_get_proc_rank_xyz, &
 s_get_proc_rank, &
 s_compute_mixlayer_thickenss, &
 s_detect_qsv, &
 s_derive_numerical_schlieren_function, &
 s_compute_speed_of_sound, &
 s_finalize_derived_variables_module

    real(wp), allocatable, dimension(:, :, :) :: gm_rho_sf !<
    !! Gradient magnitude (gm) of the density for each cell of the computational
    !! sub-domain. This variable is employed in the calculation of the numerical
    !! Schlieren function.

    !> @name Finite-difference (fd) coefficients in x-, y- and z-coordinate directions.
    !! Note that because sufficient boundary information is available for all the
    !! active coordinate directions, the centered family of the finite-difference
    !! schemes is used.
    !> @{
    real(wp), allocatable, dimension(:, :), public :: fd_coeff_x
    real(wp), allocatable, dimension(:, :), public :: fd_coeff_y
    real(wp), allocatable, dimension(:, :), public :: fd_coeff_z
    !> @}

    integer :: proc_rank_x, proc_rank_y, proc_rank_tmp, proc_rank_z

    integer, private :: flg  !<
    !! Flagging (flg) variable used to annotate the dimensionality of the dataset
    !! that is undergoing the post-process. A flag value of 1 indicates that the
    !! dataset is 3D, while a flag value of 0 indicates that it is not. This flg
    !! variable is necessary to avoid cycling through the third dimension of the
    !! flow variable(s) when the simulation is not 3D and the size of the buffer
    !! is non-zero. Note that a similar procedure does not have to be applied to
    !! the second dimension since in 1D, the buffer size is always zero.

contains

    !>  Computation of parameters, allocation procedures, and/or
        !!      any other tasks needed to properly setup the module
    impure subroutine s_initialize_derived_variables_module

        ! Allocating the gradient magnitude of the density variable provided
        ! that numerical Schlieren function is outputted during post-process
        if (schlieren_wrt) then
            allocate (gm_rho_sf(-offset_x%beg:m + offset_x%end, &
                                -offset_y%beg:n + offset_y%end, &
                                -offset_z%beg:p + offset_z%end))
        end if

        ! Allocating the variables which will store the coefficients of the
        ! centered family of finite-difference schemes. Note that sufficient
        ! space is allocated so that the coefficients up to any chosen order
        ! of accuracy may be bookkept. However, if higher than fourth-order
        ! accuracy coefficients are wanted, the formulae required to compute
        ! these coefficients will have to be implemented in the subroutine
        ! s_compute_finite_difference_coefficients.

        ! Allocating centered finite-difference coefficients in x-direction
        if (omega_wrt(2) .or. omega_wrt(3) .or. schlieren_wrt .or. liutex_wrt) then
            allocate (fd_coeff_x(-fd_number:fd_number, &
                                 -offset_x%beg:m + offset_x%end))
        end if

        ! Allocating centered finite-difference coefficients in y-direction
        if (omega_wrt(1) .or. omega_wrt(3) .or. liutex_wrt &
            .or. &
            (n > 0 .and. schlieren_wrt)) then
            allocate (fd_coeff_y(-fd_number:fd_number, &
                                 -offset_y%beg:n + offset_y%end))
        end if

        ! Allocating centered finite-difference coefficients in z-direction
        if (omega_wrt(1) .or. omega_wrt(2) .or. liutex_wrt &
            .or. &
            (p > 0 .and. schlieren_wrt)) then
            allocate (fd_coeff_z(-fd_number:fd_number, &
                                 -offset_z%beg:p + offset_z%end))
        end if

        ! Annotating the dimensionality of the dataset undergoing the post-
        ! process. A flag value of 1 indicates that the dataset is 3D, while
        ! a flag value of 0 indicates that it is not.
        if (p > 0) then
            flg = 1
        else
            flg = 0
        end if

    end subroutine s_initialize_derived_variables_module

    !>  This subroutine receives as input the specific heat ratio
        !!      function, gamma_sf, and derives from it the specific heat
        !!      ratio. The latter is stored in the derived flow quantity
        !!      storage variable, q_sf.
        !!  @param q_sf Specific heat ratio
    pure subroutine s_derive_specific_heat_ratio(q_sf)

        real(wp), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end, &
                      -offset_z%beg:p + offset_z%end), &
            intent(inout) :: q_sf

        integer :: i, j, k !< Generic loop iterators

        ! Computing specific heat ratio from specific heat ratio function
        do k = -offset_z%beg, p + offset_z%end
            do j = -offset_y%beg, n + offset_y%end
                do i = -offset_x%beg, m + offset_x%end
                    q_sf(i, j, k) = 1._wp + 1._wp/gamma_sf(i, j, k)
                end do
            end do
        end do

    end subroutine s_derive_specific_heat_ratio

    !>  This subroutine admits as inputs the specific heat ratio
        !!      function and the liquid stiffness function, gamma_sf and
        !!      pi_inf_sf, respectively. These are used to calculate the
        !!      values of the liquid stiffness, which are stored in the
        !!      derived flow quantity storage variable, q_sf.
        !!  @param q_sf Liquid stiffness
    pure subroutine s_derive_liquid_stiffness(q_sf)

        real(wp), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end, &
                      -offset_z%beg:p + offset_z%end), &
            intent(inout) :: q_sf

        integer :: i, j, k !< Generic loop iterators

        ! Calculating the values of the liquid stiffness from those of the
        ! specific heat ratio function and the liquid stiffness function
        do k = -offset_z%beg, p + offset_z%end
            do j = -offset_y%beg, n + offset_y%end
                do i = -offset_x%beg, m + offset_x%end
                    q_sf(i, j, k) = pi_inf_sf(i, j, k)/(gamma_sf(i, j, k) + 1._wp)
                end do
            end do
        end do

    end subroutine s_derive_liquid_stiffness

    !> This subroutine admits as inputs the primitive variables,
        !!      the density, the specific heat ratio function and liquid
        !!      stiffness function. It then computes from those variables
        !!      the values of the speed of sound, which are stored in the
        !!      derived flow quantity storage variable, q_sf.
        !! @param q_prim_vf Primitive variables
        !! @param q_sf Speed of sound
    pure subroutine s_derive_sound_speed(q_prim_vf, q_sf)

        type(scalar_field), &
            dimension(sys_size), &
            intent(in) :: q_prim_vf

        real(wp), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end, &
                      -offset_z%beg:p + offset_z%end), &
            intent(inout) :: q_sf

        integer :: i, j, k !< Generic loop iterators

        ! Fluid bulk modulus for alternate sound speed
        real(wp) :: blkmod1, blkmod2

        ! Computing speed of sound values from those of pressure, density,
        ! specific heat ratio function and the liquid stiffness function
        do k = -offset_z%beg, p + offset_z%end
            do j = -offset_y%beg, n + offset_y%end
                do i = -offset_x%beg, m + offset_x%end

                    ! Compute mixture sound speed
                    if (alt_soundspeed .neqv. .true.) then
                        q_sf(i, j, k) = (((gamma_sf(i, j, k) + 1._wp)* &
                                          q_prim_vf(E_idx)%sf(i, j, k) + &
                                          pi_inf_sf(i, j, k))/(gamma_sf(i, j, k)* &
                                                               rho_sf(i, j, k)))
                    else
                        blkmod1 = ((fluid_pp(1)%gamma + 1._wp)*q_prim_vf(E_idx)%sf(i, j, k) + &
                                   fluid_pp(1)%pi_inf)/fluid_pp(1)%gamma
                        blkmod2 = ((fluid_pp(2)%gamma + 1._wp)*q_prim_vf(E_idx)%sf(i, j, k) + &
                                   fluid_pp(2)%pi_inf)/fluid_pp(2)%gamma
                        q_sf(i, j, k) = (1._wp/(rho_sf(i, j, k)*(q_prim_vf(adv_idx%beg)%sf(i, j, k)/blkmod1 + &
                                                                 (1._wp - q_prim_vf(adv_idx%beg)%sf(i, j, k))/blkmod2)))
                    end if

                    if (mixture_err .and. q_sf(i, j, k) < 0._wp) then
                        q_sf(i, j, k) = 1.e-16_wp
                    else
                        q_sf(i, j, k) = sqrt(q_sf(i, j, k))
                    end if
                end do
            end do
        end do

    end subroutine s_derive_sound_speed

    !>  This subroutine derives the flux_limiter at cell boundary
        !!      i+1/2. This is an approximation because the velocity used
        !!      to determine the upwind direction is the velocity at the
        !!      cell center i instead of the contact velocity at the cell
        !!      boundary from the Riemann solver.
        !!  @param i Component indicator
        !!  @param q_prim_vf Primitive variables
        !!  @param q_sf Flux limiter
    pure subroutine s_derive_flux_limiter(i, q_prim_vf, q_sf)

        integer, intent(in) :: i

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf

        real(wp), dimension(-offset_x%beg:m + offset_x%end, &
                            -offset_y%beg:n + offset_y%end, &
                            -offset_z%beg:p + offset_z%end), &
            intent(inout) :: q_sf

        real(wp) :: top, bottom, slope !< Flux limiter calcs
        integer :: j, k, l !< Generic loop iterators

        do l = -offset_z%beg, p + offset_z%end
            do k = -offset_y%beg, n + offset_y%end
                do j = -offset_x%beg, m + offset_x%end
                    if (i == 1) then
                        if (q_prim_vf(cont_idx%end + i)%sf(j, k, l) >= 0._wp) then
                            top = q_prim_vf(adv_idx%beg)%sf(j, k, l) - &
                                  q_prim_vf(adv_idx%beg)%sf(j - 1, k, l)
                            bottom = q_prim_vf(adv_idx%beg)%sf(j + 1, k, l) - &
                                     q_prim_vf(adv_idx%beg)%sf(j, k, l)
                        else
                            top = q_prim_vf(adv_idx%beg)%sf(j + 2, k, l) - &
                                  q_prim_vf(adv_idx%beg)%sf(j + 1, k, l)
                            bottom = q_prim_vf(adv_idx%beg)%sf(j + 1, k, l) - &
                                     q_prim_vf(adv_idx%beg)%sf(j, k, l)
                        end if
                    elseif (i == 2) then
                        if (q_prim_vf(cont_idx%end + i)%sf(j, k, l) >= 0._wp) then
                            top = q_prim_vf(adv_idx%beg)%sf(j, k, l) - &
                                  q_prim_vf(adv_idx%beg)%sf(j, k - 1, l)
                            bottom = q_prim_vf(adv_idx%beg)%sf(j, k + 1, l) - &
                                     q_prim_vf(adv_idx%beg)%sf(j, k, l)
                        else
                            top = q_prim_vf(adv_idx%beg)%sf(j, k + 2, l) - &
                                  q_prim_vf(adv_idx%beg)%sf(j, k + 1, l)
                            bottom = q_prim_vf(adv_idx%beg)%sf(j, k + 1, l) - &
                                     q_prim_vf(adv_idx%beg)%sf(j, k, l)
                        end if
                    else
                        if (q_prim_vf(cont_idx%end + i)%sf(j, k, l) >= 0._wp) then
                            top = q_prim_vf(adv_idx%beg)%sf(j, k, l) - &
                                  q_prim_vf(adv_idx%beg)%sf(j, k, l - 1)
                            bottom = q_prim_vf(adv_idx%beg)%sf(j, k, l + 1) - &
                                     q_prim_vf(adv_idx%beg)%sf(j, k, l)
                        else
                            top = q_prim_vf(adv_idx%beg)%sf(j, k, l + 2) - &
                                  q_prim_vf(adv_idx%beg)%sf(j, k, l + 1)
                            bottom = q_prim_vf(adv_idx%beg)%sf(j, k, l + 1) - &
                                     q_prim_vf(adv_idx%beg)%sf(j, k, l)
                        end if
                    end if

                    if (abs(top) < 1.e-8_wp) top = 0._wp
                    if (abs(bottom) < 1.e-8_wp) bottom = 0._wp

                    if (f_approx_equal(top, bottom)) then
                        slope = 1._wp
                        !       ELSEIF((top == 0._wp .AND. bottom /= 0._wp) &
                        !               .OR.            &
                        !           (bottom == 0._wp .AND. top /= 0._wp)) THEN
                        !           slope = 0._wp
                    else
                        slope = (top*bottom)/(bottom**2._wp + 1.e-16_wp)
                    end if

                    ! Flux limiter function
                    if (flux_lim == 1) then ! MINMOD (MM)
                        q_sf(j, k, l) = max(0._wp, min(1._wp, slope))
                    elseif (flux_lim == 2) then ! MUSCL (MC)
                        q_sf(j, k, l) = max(0._wp, min(2._wp*slope, 5.e-1_wp*(1._wp + slope), 2._wp))
                    elseif (flux_lim == 3) then ! OSPRE (OP)
                        q_sf(j, k, l) = (15.e-1_wp*(slope**2._wp + slope))/(slope**2._wp + slope + 1._wp)
                    elseif (flux_lim == 4) then ! SUPERBEE (SB)
                        q_sf(j, k, l) = max(0._wp, min(1._wp, 2._wp*slope), min(slope, 2._wp))
                    elseif (flux_lim == 5) then ! SWEBY (SW) (beta = 1.5)
                        q_sf(j, k, l) = max(0._wp, min(15.e-1_wp*slope, 1._wp), min(slope, 15.e-1_wp))
                    elseif (flux_lim == 6) then ! VAN ALBADA (VA)
                        q_sf(j, k, l) = (slope**2._wp + slope)/(slope**2._wp + 1._wp)
                    elseif (flux_lim == 7) then ! VAN LEER (VL)
                        q_sf(j, k, l) = (abs(slope) + slope)/(1._wp + abs(slope))
                    end if
                end do
            end do
        end do
    end subroutine s_derive_flux_limiter

    !>  Computes the solution to the linear system Ax=b w/ sol = x
        !!  @param A Input matrix
        !!  @param b right-hane-side
        !!  @param sol Solution
        !!  @param ndim Problem size
    pure subroutine s_solve_linear_system(A, b, sol, ndim)

        integer, intent(in) :: ndim
        real(wp), dimension(ndim, ndim), intent(inout) :: A
        real(wp), dimension(ndim), intent(inout) :: b
        real(wp), dimension(ndim), intent(out) :: sol

        !EXTERNAL DGESV

        integer :: i, j, k

        ! Solve linear system using own linear solver (Thomson/Darter/Comet/Stampede)
        ! Forward elimination
        do i = 1, ndim
            ! Pivoting
            j = i - 1 + maxloc(abs(A(i:ndim, i)), 1)
            sol = A(i, :)
            A(i, :) = A(j, :)
            A(j, :) = sol
            sol(1) = b(i)
            b(i) = b(j)
            b(j) = sol(1)
            ! Elimination
            b(i) = b(i)/A(i, i)
            A(i, :) = A(i, :)/A(i, i)
            do k = i + 1, ndim
                b(k) = b(k) - A(k, i)*b(i)
                A(k, :) = A(k, :) - A(k, i)*A(i, :)
            end do
        end do

        ! Backward substitution
        do i = ndim, 1, -1
            sol(i) = b(i)
            do k = i + 1, ndim
                sol(i) = sol(i) - A(i, k)*sol(k)
            end do
        end do

    end subroutine s_solve_linear_system

    !>  This subroutine receives as inputs the indicator of the
        !!      component of the vorticity that should be outputted and
        !!      the primitive variables. From those inputs, it proceeds
        !!      to calculate values of the desired vorticity component,
        !!      which are subsequently stored in derived flow quantity
        !!      storage variable, q_sf.
        !!  @param i Vorticity component indicator
        !!  @param q_prim_vf Primitive variables
        !!  @param q_sf Vorticity component
    pure subroutine s_derive_vorticity_component(i, q_prim_vf, q_sf)

        integer, intent(in) :: i

        type(scalar_field), &
            dimension(sys_size), &
            intent(in) :: q_prim_vf

        real(wp), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end, &
                      -offset_z%beg:p + offset_z%end), &
            intent(inout) :: q_sf

        integer :: j, k, l, r !< Generic loop iterators

        ! Computing the vorticity component in the x-coordinate direction
        if (i == 1) then
            do l = -offset_z%beg, p + offset_z%end
                do k = -offset_y%beg, n + offset_y%end
                    do j = -offset_x%beg, m + offset_x%end

                        q_sf(j, k, l) = 0._wp

                        do r = -fd_number, fd_number
                            if (grid_geometry == 3) then
                                q_sf(j, k, l) = &
                                    q_sf(j, k, l) + 1._wp/y_cc(k)* &
                                    (fd_coeff_y(r, k)*y_cc(r + k)* &
                                     q_prim_vf(mom_idx%end)%sf(j, r + k, l) &
                                     - fd_coeff_z(r, l)* &
                                     q_prim_vf(mom_idx%beg + 1)%sf(j, k, r + l))
                            else
                                q_sf(j, k, l) = &
                                    q_sf(j, k, l) + fd_coeff_y(r, k)* &
                                    q_prim_vf(mom_idx%end)%sf(j, r + k, l) &
                                    - fd_coeff_z(r, l)* &
                                    q_prim_vf(mom_idx%beg + 1)%sf(j, k, r + l)
                            end if
                        end do

                    end do
                end do
            end do

            ! Computing the vorticity component in the y-coordinate direction
        elseif (i == 2) then
            do l = -offset_z%beg, p + offset_z%end
                do k = -offset_y%beg, n + offset_y%end
                    do j = -offset_x%beg, m + offset_x%end

                        q_sf(j, k, l) = 0._wp

                        do r = -fd_number, fd_number
                            if (grid_geometry == 3) then
                                q_sf(j, k, l) = &
                                    q_sf(j, k, l) + fd_coeff_z(r, l)/y_cc(k)* &
                                    q_prim_vf(mom_idx%beg)%sf(j, k, r + l) &
                                    - fd_coeff_x(r, j)* &
                                    q_prim_vf(mom_idx%end)%sf(r + j, k, l)
                            else
                                q_sf(j, k, l) = &
                                    q_sf(j, k, l) + fd_coeff_z(r, l)* &
                                    q_prim_vf(mom_idx%beg)%sf(j, k, r + l) &
                                    - fd_coeff_x(r, j)* &
                                    q_prim_vf(mom_idx%end)%sf(r + j, k, l)
                            end if
                        end do

                    end do
                end do
            end do

            ! Computing the vorticity component in the z-coordinate direction
        else
            do l = -offset_z%beg, p + offset_z%end
                do k = -offset_y%beg, n + offset_y%end
                    do j = -offset_x%beg, m + offset_x%end

                        q_sf(j, k, l) = 0._wp

                        do r = -fd_number, fd_number
                            q_sf(j, k, l) = &
                                q_sf(j, k, l) + fd_coeff_x(r, j)* &
                                q_prim_vf(mom_idx%beg + 1)%sf(r + j, k, l) &
                                - fd_coeff_y(r, k)* &
                                q_prim_vf(mom_idx%beg)%sf(j, r + k, l)
                        end do

                    end do
                end do
            end do
        end if

    end subroutine s_derive_vorticity_component

    !> This subroutine gets as inputs the primitive variables. From those
        !!      inputs, it proceeds to calculate the value of the Q_M
        !!      function, which are subsequently stored in the derived flow
        !!      quantity storage variable, q_sf.
        !!  @param q_prim_vf Primitive variables
        !!  @param q_sf Q_M
    pure subroutine s_derive_qm(q_prim_vf, q_sf)
        type(scalar_field), &
            dimension(sys_size), &
            intent(in) :: q_prim_vf

        real(wp), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end, &
                      -offset_z%beg:p + offset_z%end), &
            intent(inout) :: q_sf

        real(wp), &
            dimension(1:3, 1:3) :: q_jacobian_sf, S, S2, O, O2

        real(wp) :: trS, Q, IIS
        integer :: j, k, l, r, jj, kk !< Generic loop iterators

        do l = -offset_z%beg, p + offset_z%end
            do k = -offset_y%beg, n + offset_y%end
                do j = -offset_x%beg, m + offset_x%end

                    ! Get velocity gradient tensor
                    q_jacobian_sf(:, :) = 0._wp

                    do r = -fd_number, fd_number
                        do jj = 1, 3
                            ! d()/dx
                            q_jacobian_sf(jj, 1) = &
                                q_jacobian_sf(jj, 1) + &
                                fd_coeff_x(r, j)* &
                                q_prim_vf(mom_idx%beg + jj - 1)%sf(r + j, k, l)
                            ! d()/dy
                            q_jacobian_sf(jj, 2) = &
                                q_jacobian_sf(jj, 2) + &
                                fd_coeff_y(r, k)* &
                                q_prim_vf(mom_idx%beg + jj - 1)%sf(j, r + k, l)
                            ! d()/dz
                            q_jacobian_sf(jj, 3) = &
                                q_jacobian_sf(jj, 3) + &
                                fd_coeff_z(r, l)* &
                                q_prim_vf(mom_idx%beg + jj - 1)%sf(j, k, r + l)
                        end do
                    end do

                    ! Decompose J into asymmetric matrix, S, and a skew-symmetric matrix, O
                    do jj = 1, 3
                        do kk = 1, 3
                            S(jj, kk) = 0.5_wp* &
                                        (q_jacobian_sf(jj, kk) + q_jacobian_sf(kk, jj))
                            O(jj, kk) = 0.5_wp* &
                                        (q_jacobian_sf(jj, kk) - q_jacobian_sf(kk, jj))
                        end do
                    end do

                    ! Compute S2 = S*S'
                    do jj = 1, 3
                        do kk = 1, 3
                            O2(jj, kk) = O(jj, 1)*O(kk, 1) + &
                                         O(jj, 2)*O(kk, 2) + &
                                         O(jj, 3)*O(kk, 3)
                            S2(jj, kk) = S(jj, 1)*S(kk, 1) + &
                                         S(jj, 2)*S(kk, 2) + &
                                         S(jj, 3)*S(kk, 3)
                        end do
                    end do

                    ! Compute Q
                    Q = 0.5_wp*((O2(1, 1) + O2(2, 2) + O2(3, 3)) - &
                                (S2(1, 1) + S2(2, 2) + S2(3, 3)))
                    trS = S(1, 1) + S(2, 2) + S(3, 3)
                    IIS = 0.5_wp*((S(1, 1) + S(2, 2) + S(3, 3))**2 - &
                                  (S2(1, 1) + S2(2, 2) + S2(3, 3)))
                    q_sf(j, k, l) = Q + IIS

                end do
            end do
        end do

    end subroutine s_derive_qm

    !> This subroutine gets as inputs the primitive variables. From those
        !!      inputs, it proceeds to calculate the Liutex vector and its
        !!      magnitude based on Xu et al. (2019).
        !!  @param q_prim_vf Primitive variables
        !!  @param liutex_mag Liutex magnitude
        !!  @param liutex_axis Liutex axis
    impure subroutine s_derive_liutex(q_prim_vf, y_idx_beg, y_idx_end, &
                                    liutex_mag, liutex_axis, liutex_rrs, liutex_core, &
                                    omega, omega_axis, omega_perp, &
                                    vort_stretch, vort_stretch_proj, vort_stretch_res, &
                                    A_rr, A_ps, A_ns, A_sr)
        integer, parameter :: ndim = 3
        type(scalar_field), &
            dimension(sys_size), &
            intent(in) :: q_prim_vf

        integer, intent(in) :: y_idx_beg, y_idx_end

        real(wp), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end, &
                      -offset_z%beg:p + offset_z%end), &
            intent(out) :: liutex_mag, liutex_rrs, liutex_core, omega_axis, omega_perp

        real(wp), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end, &
                      -offset_z%beg:p + offset_z%end, ndim), &
            intent(out) :: liutex_axis !< Liutex rigid rotation axis

        real(wp), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end, &
                      -offset_z%beg:p + offset_z%end, ndim), &
            intent(out) :: omega, vort_stretch !< Vorticity

        real(wp), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end, &
                      -offset_z%beg:p + offset_z%end), &
            intent(out) :: vort_stretch_proj, vort_stretch_res, A_rr, A_ps, A_ns, A_sr

        real(wp), dimension(ndim, ndim) :: A_S, A_W
        real(wp), dimension(ndim) :: vort_ps
        real(wp) :: S2, W2

        character, parameter :: ivl = 'N' !< compute left eigenvectors
        character, parameter :: ivr = 'V' !< compute right eigenvectors
        real(wp), dimension(ndim, ndim) :: vgt !< velocity gradient tensor
        real(wp), dimension(ndim) :: lr, li !< real and imaginary parts of eigenvalues
        real(wp), dimension(ndim, ndim) :: vl, vr !< left and right eigenvectors
        integer, parameter :: lwork = 4*ndim !< size of work array (4*ndim recommended)
        real(wp), dimension(lwork) :: work !< work array
        integer :: info

        real(wp), dimension(ndim) :: eigvec !< real eigenvector
        real(wp) :: eigvec_mag !< magnitude of real eigenvector
        real(wp) :: omega_proj !< projection of vorticity on real eigenvector
        real(wp) :: lrr, lcr, lci !< imaginary part of complex eigenvalue
        real(wp) :: alpha, beta

        real(wp), dimension(3) :: grad_liutex_mag
        real(wp) :: liutex_core_val

        real(wp), dimension(3) :: tmp

        integer :: j, k, l, r, i !< Generic loop iterators
        integer :: idx

        omega = 0._wp
        liutex_mag = 0._wp
        liutex_axis = 0._wp
        liutex_rrs = 0._wp 
        vort_stretch_proj = 0._wp 
        vort_stretch_res = 0._wp 
        A_rr = 0._wp
        A_ps = 0._wp
        A_ns = 0._wp
        A_sr = 0._wp
        do k = -offset_y%beg, n + offset_y%end 
            if (y_cc(k) >= y_cc_glb(y_idx_beg) .and. y_cc(k) <= y_cc_glb(y_idx_end)) then
                do l = -offset_z%beg, p + offset_z%end
                    do j = -offset_x%beg, m + offset_x%end

                        ! Get velocity gradient tensor (VGT)
                        vgt(:, :) = 0._wp

                        do r = -fd_number, fd_number
                            do i = 1, 3
                                ! d()/dx
                                vgt(i, 1) = &
                                    vgt(i, 1) + &
                                    fd_coeff_x(r, j)* &
                                    q_prim_vf(mom_idx%beg + i - 1)%sf(r + j, k, l)
                                ! d()/dy
                                vgt(i, 2) = &
                                    vgt(i, 2) + &
                                    fd_coeff_y(r, k)* &
                                    q_prim_vf(mom_idx%beg + i - 1)%sf(j, r + k, l)
                                ! d()/dz
                                vgt(i, 3) = &
                                    vgt(i, 3) + &
                                    fd_coeff_z(r, l)* &
                                    q_prim_vf(mom_idx%beg + i - 1)%sf(j, k, r + l)
                            end do
                        end do

                        ! Compute vorticity
                        omega(j, k, l, 1) = vgt(3,2) - vgt(2,3)
                        omega(j, k, l, 2) = vgt(1,3) - vgt(3,1)
                        omega(j, k, l, 3) = vgt(2,1) - vgt(1,2)

                        ! Compute vortex stretching term
                        do r = 1, 3
                            vort_stretch(j, k, l, r) = omega(j, k, l, 1)*vgt(r,1) &
                                                    + omega(j, k, l, 2)*vgt(r,2) &
                                                    + omega(j, k, l, 3)*vgt(r,3)
                        end do

                        
                        A_S = 0.5_wp*(vgt + transpose(vgt))
                        A_W = 0.5_wp*(vgt - transpose(vgt))
                        S2 = 0._wp
                        W2 = 0._wp
                        do i = 1, 3
                            do r = 1, 3
                                S2 = S2 + A_S(i, r)**2._wp
                                W2 = W2 + A_W(i, r)**2._wp
                            end do
                        end do

                        ! Call appropriate LAPACK routine based on precision
#ifdef MFC_SINGLE_PRECISION
                        call sgeev(ivl, ivr, ndim, vgt, ndim, lr, li, vl, ndim, vr, ndim, work, lwork, info)
#else
                        call dgeev(ivl, ivr, ndim, vgt, ndim, lr, li, vl, ndim, vr, ndim, work, lwork, info)
#endif

                        ! Find real eigenvector
                        idx = 1
                        do r = 2, 3
                            if (abs(li(r)) < abs(li(idx))) then
                                idx = r
                            end if
                        end do
                        eigvec = vr(:, idx)

                        ! Normalize real eigenvector if it is effectively non-zero
                        eigvec_mag = sqrt(eigvec(1)**2._wp &
                                          + eigvec(2)**2._wp &
                                          + eigvec(3)**2._wp)
                        if (eigvec_mag > sgm_eps) then
                            eigvec = eigvec/eigvec_mag
                        else
                            eigvec = 0._wp
                        end if

                        ! Compute vorticity projected on the eigenvector
                        omega_proj = sum(omega(j, k, l, :)*eigvec(:))
                        omega_axis(j, k, l) = omega_proj

                        ! Compute vorticity perpendicular to the eigenvector
                        tmp(:) = omega(j, k, l, :) - omega_proj*eigvec(:)
                        omega_perp(j, k, l) = sqrt(sum(tmp**2._wp))

                        ! As eigenvector can have +/- signs, we can choose the sign
                        ! so that omega_proj is positive
                        if (omega_proj < 0._wp) then
                            eigvec = -eigvec
                            omega_proj = -omega_proj
                        end if

                        ! Find real and imaginary part of complex eigenvalue
                        lrr = lr(idx)
                        lcr = lr(mod(idx, 3) + 1)
                        lci = li(mod(idx, 3) + 1)
                        
                        ! Compute Liutex magnitude
                        alpha = 0.25_wp*omega_proj**2._wp - lci**2._wp
                        if (alpha > 0._wp) then
                            alpha = sqrt(alpha)
                        else
                            alpha = 0._wp
                        end if
                        beta = 0.5_wp*omega_proj

                        ! Compute Liutex magnitude
                        liutex_mag(j, k, l) = 2._wp*(beta - alpha)
                        ! Compute relative rotation strength
                        liutex_rrs(j, k, l) = beta**2._wp / (beta**2._wp + alpha**2._wp + lcr**2._wp + 0.5_wp*lrr**2._wp + sgm_eps)
                        ! Compute Liutex axis
                        liutex_axis(j, k, l, :) = eigvec(:)

                        ! Compute projection of vortex stretching term on Liutex axis
                        vort_stretch_proj(j, k, l) = &
                            abs(liutex_axis(j, k, l, 1)*vort_stretch(j, k, l, 1) + &
                                liutex_axis(j, k, l, 2)*vort_stretch(j, k, l, 2) + &
                                liutex_axis(j, k, l, 3)*vort_stretch(j, k, l, 3))

                        vort_stretch_res(j, k, l) = &
                            sqrt((vort_stretch(j, k, l, 1) - liutex_axis(j, k, l, 1)*vort_stretch_proj(j, k, l))**2._wp + &
                                (vort_stretch(j, k, l, 2) - liutex_axis(j, k, l, 2)*vort_stretch_proj(j, k, l))**2._wp + &
                                (vort_stretch(j, k, l, 3) - liutex_axis(j, k, l, 3)*vort_stretch_proj(j, k, l))**2._wp)

                        ! Strength
                        vort_ps(1) = omega(j, k, l, 1) - liutex_mag(j, k, l)*liutex_axis(j, k, l, 1)
                        vort_ps(2) = omega(j, k, l, 2) - liutex_mag(j, k, l)*liutex_axis(j, k, l, 2)
                        vort_ps(3) = omega(j, k, l, 3) - liutex_mag(j, k, l)*liutex_axis(j, k, l, 3)
                        A_rr(j, k, l) = 0.5_wp * liutex_mag(j, k, l)**2._wp
                        A_ps(j, k, l) = vort_ps(1)**2._wp + vort_ps(2)**2._wp + vort_ps(3)**2._wp
                        A_ns(j, k, l) = S2 - 0.5_wp*A_ps(j, k, l)
                        A_sr(j, k, l) = W2 - 0.5_wp*A_ps(j, k, l) - 0.5_wp*A_rr(j, k, l)
                    end do
                end do
            end if
        end do

    end subroutine s_derive_liutex

    impure subroutine s_compute_mixlayer_thickenss(vel1_loc, mixlayer_thickness, mixlayer_idx_beg, mixlayer_idx_end)
        real(wp), intent(in), dimension(0:m, 0:n, 0:p) :: vel1_loc
        real(wp), intent(out) :: mixlayer_thickness
        integer, intent(out) :: mixlayer_idx_beg, mixlayer_idx_end
        real(wp), dimension(0:n) :: vel1_loc_sum
        real(wp), dimension(0:n_glb) :: vel1_sum, vel1_glb_avg
        real(wp), dimension(0:n_glb) :: y_loc
        integer :: j, jj, ierr 

        ! Compute local average field for each processor
        vel1_loc_sum = sum(sum(vel1_loc, 3), 1)

        ! Compute global average field
        vel1_sum = 0._wp
        do j = 0, n
            jj = j + proc_rank_y*(n + 1)
            vel1_sum(jj) = vel1_sum(jj) + vel1_loc_sum(j)
        end do
        call MPI_ALLREDUCE(vel1_sum, vel1_glb_avg(0:n_glb), n_glb + 1, mpi_p, &
                           MPI_SUM, MPI_COMM_WORLD, ierr)
        vel1_glb_avg = vel1_glb_avg / ((m_glb + 1) * (p_glb + 1))

        ! Compute mixlayer_thickness
        if (proc_rank == 0) then 
            y_loc = y_cc_glb(0:n_glb)
            do j = 0, n_glb
                if (vel1_glb_avg(j) < -0.99_wp) y_loc(j) = 0._wp
                if (vel1_glb_avg(j) >  0.99_wp) y_loc(j) = 0._wp
            end do
            mixlayer_idx_beg = minloc(y_loc, DIM=1)
            mixlayer_idx_end = maxloc(y_loc, DIM=1)
            mixlayer_thickness = y_loc(mixlayer_idx_end) - y_loc(mixlayer_idx_beg)
        end if
        call MPI_BCAST(mixlayer_thickness, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(mixlayer_idx_beg, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(mixlayer_idx_end, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
    end subroutine s_compute_mixlayer_thickenss

    impure subroutine s_detect_qsv(liutex_mag, liutex_axis, omega_axis, omega_perp, A_rr, A_ps, y_idx_beg, y_idx_end, qsv_info, q_sf_group)
        real(wp), dimension(-offset_x%beg:m + offset_x%end, &
                            -offset_y%beg:n + offset_y%end, &
                            -offset_z%beg:p + offset_z%end), intent(in) :: liutex_mag, omega_axis, omega_perp, A_rr, A_ps
        real(wp), dimension(-offset_x%beg:m + offset_x%end, &
                            -offset_y%beg:n + offset_y%end, &
                            -offset_z%beg:p + offset_z%end, 3), intent(in) :: liutex_axis

        integer, intent(in) :: y_idx_beg, y_idx_end

        real(wp), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end, &
                      -offset_z%beg:p + offset_z%end, 5), intent(out) :: qsv_info

        real(wp), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end, &
                      -offset_z%beg:p + offset_z%end), intent(inout) :: q_sf_group
                      
        logical, dimension(0:m, 0:n, 0:p, 5) :: qsv_flag
        integer, dimension(0:m, 0:n, 0:p) :: qsv_group
        integer :: qsv_flag_count
        logical, dimension(-1:m + 1, -1:n + 1, -1:p + 1) :: qsv_flag_padded
        real(wp) :: theta1, theta2
        integer :: i, j, k, l, ierr

        ! Initialization
        qsv_flag = .false.
        qsv_flag_count = 1
        qsv_info = 0._wp
        q_sf_group = 0._wp

        ! (1,2,3) Initial filtering
        if (proc_rank == 0) print *, "criteria 1-3"
        do j = 0, n
            if (y_cc(j) >= y_cc_glb(y_idx_beg) .and. y_cc(j) <= y_cc_glb(y_idx_end)) then
                do k = 0, p
                    do i = 0, m
                        ! (1) liutex_mag
                        if (liutex_mag(i, j, k) > 1.e-10_wp) then
                            qsv_flag(i, j, k, 1) = .true.
                            qsv_info(i, j, k, 1) = 1._wp
                        end if

                        ! (2) liutex axis
                        if (qsv_flag(i, j, k, 1)) then
                            theta1 = atan(liutex_axis(i, j, k, 2) / liutex_axis(i, j, k, 1)) / pi * 180._wp
                            theta2 = atan(liutex_axis(i, j, k, 3) / liutex_axis(i, j, k, 1)) / pi * 180._wp
                            if (theta1 > 0._wp .and. theta1 < 90._wp .and. &
                                theta2 > -45._wp .and. theta2 < 45._wp) then
                                qsv_flag(i, j, k, 2) = .true.
                                qsv_info(i, j, k, 2) = 1._wp
                            end if
                        end if

                        ! (3) A_rr > A_ps
                        if (qsv_flag(i, j, k, 2) .and. &
                            A_rr(i, j, k) > A_ps(i, j, k)) then
                            qsv_flag(i, j, k, 3) = .true.
                            qsv_info(i, j, k, 3) = 1._wp
                        end if

                    end do
                end do
            end if
        end do

        ! (4) Remove isolated single point in each yz plane
        if (proc_rank == 0) print *, "criteria 4"
        do j = 0, n
            if (y_cc(j) >= y_cc_glb(y_idx_beg) .and. y_cc(j) <= y_cc_glb(y_idx_end)) then
                do k = 0, p
                    do i = 0, m
                        if (qsv_flag(i, j, k, 3)) then 
                          call s_remove_isolated_single_point(i, j, k, qsv_flag(:, :, :, 3), qsv_flag(:, :, :, 4))
                          if (qsv_flag(i, j, k, 4)) qsv_info(i, j, k, 4) = 1._wp
                        end if
                    end do
                end do
            end if
        end do

        if (proc_rank == 0) print *, "criteria 5"
        call s_classify_groups(qsv_flag, qsv_info(0:m, 0:n, 0:p, 5), qsv_group, y_idx_beg, y_idx_end)
        q_sf_group(0:m, 0:n, 0:p) = real(qsv_group, wp)
        if (proc_rank == 0) print *, "s_detect_qsv done"

    end subroutine s_detect_qsv

    impure subroutine s_classify_groups(qsv_flag, qsv_info, qsv_group, y_idx_beg, y_idx_end)
        logical, dimension(0:m, 0:n, 0:p, 5), intent(inout) :: qsv_flag
        real(wp), dimension(0:m, 0:n, 0:p), intent(out) :: qsv_info
        integer, dimension(0:m, 0:n, 0:p), intent(out) :: qsv_group
        integer, intent(in) :: y_idx_beg, y_idx_end
        integer :: id_qsv_group, id_qsv_group_max

        logical, dimension(:), allocatable :: id_qsv_group_mask, id_qsv_group_mask_glb
        integer, dimension(:), allocatable :: num_qsv_group_member, num_qsv_group_member_glb
        logical, dimension(0:m, 0:n, 0:p) :: qsv_group_mask, qsv_merge_x, qsv_merge_z
        real(wp), dimension(0:m, 0:n, 0:p) :: coord_x, coord_y, coord_z
        real(wp), dimension(0:m, 0:n, 0:p) :: x_mask, y_mask, z_mask
        real(wp), dimension(0:m, 0:n, 0:p) :: x_centered, y_centered, z_centered
        real(wp) :: x_mean, y_mean, z_mean
        real(wp) :: x_mean_glb, y_mean_glb, z_mean_glb

        integer, parameter :: ndim = 3
        character, parameter :: jobz = 'V' !< compute eigenvectors
        character, parameter :: uplo = 'U' !< Upper triangular matrix is stored
        real(wp), dimension(ndim) :: eigval !< eigenvalues
        real(wp), dimension(ndim, ndim) :: eigvec !< eigenvectors
        integer, parameter :: lwork = 4*ndim !< size of work array (4*ndim recommended)
        real(wp), dimension(lwork) :: work !< work array
        integer :: info
        real(wp), dimension(ndim) :: mean
        real(wp), dimension(ndim, ndim) :: cov, cov_glb
        real(wp), dimension(ndim) :: pca_axis
        real(wp) :: theta1, theta2
        
        integer :: ierr
        integer :: i, j, k, l

        ! Grid
        coord_x = spread(spread(x_cc(0:m), dim=2, ncopies=n + 1), dim=3, ncopies=p + 1)
        coord_y = spread(spread(y_cc(0:n), dim=1, ncopies=m + 1), dim=3, ncopies=p + 1)
        coord_z = spread(spread(z_cc(0:p), dim=1, ncopies=m + 1), dim=2, ncopies=n + 1)

        ! Grouping connected points
        if (proc_rank == 0) print *, "grouping"
        call s_identify_connected_groups(qsv_flag, qsv_group, qsv_merge_x, qsv_merge_z, id_qsv_group, y_idx_beg, y_idx_end)
        ! Find max group id
        call MPI_ALLREDUCE(id_qsv_group, id_qsv_group_max, 1, mpi_integer, MPI_MAX, MPI_COMM_WORLD, ierr)
        allocate (id_qsv_group_mask(0:id_qsv_group_max))
        allocate (id_qsv_group_mask_glb(0:id_qsv_group_max))
        id_qsv_group_mask = .false.
        do l = proc_rank + num_procs, id_qsv_group, num_procs
          id_qsv_group_mask(l) = .true.
        end do
        call MPI_ALLREDUCE(id_qsv_group_mask, id_qsv_group_mask_glb, id_qsv_group_max, mpi_logical, MPI_LOR, MPI_COMM_WORLD, ierr)
        if (proc_rank == 0) print *, "grouping done"
    
        if (proc_rank == 0) print *, "counting"
        allocate (num_qsv_group_member(0:id_qsv_group_max))
        allocate (num_qsv_group_member_glb(0:id_qsv_group_max))
        num_qsv_group_member = 0
        do l = num_procs, id_qsv_group_max
          if (id_qsv_group_mask_glb(l)) num_qsv_group_member(l) = count(qsv_group == l)
        end do
        call MPI_ALLREDUCE(num_qsv_group_member, num_qsv_group_member_glb, id_qsv_group_max, mpi_integer, MPI_SUM, MPI_COMM_WORLD, ierr)
        if (proc_rank == 0) print *, "counting done", id_qsv_group_max

        ! Perform PCA
        if (proc_rank == 0) print *, "PCA"
        do l = num_procs, id_qsv_group_max
          if (id_qsv_group_mask_glb(l) .and. num_qsv_group_member_glb(l) > 3) then
            if (proc_rank == 0) print *, "group", l, num_qsv_group_member_glb(l)
            ! Mask
            qsv_group_mask = (qsv_group == l)
            
            !
            x_mask = 0._wp
            y_mask = 0._wp
            z_mask = 0._wp
            where (qsv_group_mask)
              x_mask = coord_x
              y_mask = coord_y
              z_mask = coord_z
            end where
            where (qsv_group_mask .and. qsv_merge_x) x_mask = x_mask + (x_cb_glb(m_glb) - x_cb_glb(-1))
            where (qsv_group_mask .and. qsv_merge_z) z_mask = z_mask + (z_cb_glb(p_glb) - z_cb_glb(-1))

            ! if (l == 286) then
            !   do k = 0, p
            !     do j = 0, n
            !       do i = 0, m
            !         write(100+proc_rank,*) proc_rank_x, proc_rank_y, proc_rank_z, i, j, k, qsv_group(i, j, k), qsv_group_mask(i, j, k), qsv_merge_x(i, j, k), qsv_merge_z(i, j, k), x_mask(i, j, k), y_mask(i, j, k), z_mask(i, j, k), coord_x(i, j, k), coord_y(i, j, k), coord_z(i, j, k)
            !       end do
            !     end do
            !   end do
            ! end if

            ! Compute mean
            x_mean = sum(x_mask)
            y_mean = sum(y_mask)
            z_mean = sum(z_mask)
            call MPI_ALLREDUCE(x_mean, x_mean_glb, 1, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
            call MPI_ALLREDUCE(y_mean, y_mean_glb, 1, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
            call MPI_ALLREDUCE(z_mean, z_mean_glb, 1, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
            x_mean_glb = x_mean_glb / real(num_qsv_group_member_glb(l), wp)
            y_mean_glb = y_mean_glb / real(num_qsv_group_member_glb(l), wp)
            z_mean_glb = z_mean_glb / real(num_qsv_group_member_glb(l), wp)


            ! Center the data by subtracting mean
            x_centered = 0._wp
            y_centered = 0._wp
            z_centered = 0._wp
            where(qsv_group_mask)
              x_centered = x_mask - x_mean_glb
              y_centered = y_mask - y_mean_glb
              z_centered = z_mask - z_mean_glb
            end where

            ! Compute covariance matrix
            cov(1, 1) = sum(x_centered*x_centered)
            cov(1, 2) = sum(x_centered*y_centered)
            cov(1, 3) = sum(x_centered*z_centered)
            cov(2, 2) = sum(y_centered*y_centered)
            cov(2, 3) = sum(y_centered*z_centered)
            cov(3, 3) = sum(z_centered*z_centered)
            call MPI_ALLREDUCE(cov, cov_glb, ndim*ndim, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
            cov_glb = cov_glb / real(num_qsv_group_member_glb(l), wp)

            ! Compute eigenvalues and eigenvectors of covariance matrix
#ifdef MFC_SINGLE_PRECISION
            call ssyev(jobz, uplo, ndim, cov_glb, ndim, eigval, work, lwork, info)
#else
            call dsyev(jobz, uplo, ndim, cov_glb, ndim, eigval, work, lwork, info)
#endif

            ! Most significant eigenvector
            pca_axis = cov_glb(:, ndim)

            !
            theta1 = atan(pca_axis(2) / pca_axis(1)) / pi * 180._wp
            theta2 = atan(pca_axis(3) / pca_axis(1)) / pi * 180._wp
            ! if (l == 286) g(99, *) proc_rank_x, proc_rank_z, l, num_qsv_group_member_glb(l), eigval(ndim), pca_axis, theta1, theta2
            if (theta1 > 0._wp .and. theta1 < 90._wp .and. &
                theta2 > -45._wp .and. theta2 < 45._wp) then
                where (qsv_group_mask) qsv_flag(:, :, :, 5) = .true.
                where (qsv_group_mask) qsv_info = 1._wp
            end if
          end if
        end do

    end subroutine s_classify_groups
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    impure subroutine s_identify_connected_groups(qsv_flag, qsv_group, qsv_merge_x, qsv_merge_z, id_qsv_group, y_idx_beg, y_idx_end)
        logical, dimension(0:m, 0:n, 0:p, 5), intent(inout) :: qsv_flag
        integer, dimension(0:m, 0:n, 0:p), intent(out) :: qsv_group
        logical, dimension(0:m, 0:n, 0:p), intent(out) :: qsv_merge_x, qsv_merge_z
        integer, intent(out) :: id_qsv_group
        integer, intent(in) :: y_idx_beg, y_idx_end

        integer, dimension(-1:m + 1, -1:n + 1, -1:p + 1) :: qsv_group_padded
        logical, dimension(-1:m + 1, -1:n + 1, -1:p + 1) :: qsv_merge_x_padded
        integer :: id_qsv_group_init
        integer :: ibase, jbase, kbase, i, j, k, l, ii, jj, kk
        integer :: pr_neg, pr_pos
        integer :: ierr
        logical :: breakpoint

        id_qsv_group_init = proc_rank
        id_qsv_group = id_qsv_group_init
        qsv_group = id_qsv_group_init
        qsv_group_padded = id_qsv_group_init

        ! Find connected regions inside the domain
        if (proc_rank == 0) print *, "Find connected regions inside the domain"
        do jbase = 0, n
          if (y_cc(jbase) >= y_cc_glb(y_idx_beg) .and. y_cc(jbase) <= y_cc_glb(y_idx_end)) then
            do kbase = 0, p
              do ibase = 0, m
                if (qsv_flag(ibase, jbase, kbase, 4)) then
                  ! Check surrounding points
                  do k = kbase - 1, kbase + 1
                    do j = jbase - 1, jbase + 1
                      do i = ibase - 1, ibase + 1
                        if (i >= 0 .and. i <= m .and. &
                            j >= 0 .and. j <= n .and. &
                            k >= 0 .and. k <= p) then
                          ! if the connected point is ...
                          if (.not. (i == ibase .and. j == jbase .and. k == kbase) & ! not the base point,
                              .and. qsv_flag(i, j, k, 4)) then                       ! a qsv candidate

                              ! if both (ibase, jbase, kbase) and (i, j, k) are not grouped yet
                              if (qsv_group(i, j, k) == id_qsv_group_init .and. &
                                  qsv_group(ibase, jbase, kbase) == id_qsv_group_init) then
                                  id_qsv_group = id_qsv_group + num_procs
                                  qsv_group(i, j, k) = id_qsv_group
                                  qsv_group(ibase, jbase, kbase) = id_qsv_group

                              ! if both (ibase, jbase, kbase) and (i, j, k) are grouped, but in different groups
                              else if (qsv_group(i, j, k) /= id_qsv_group_init .and. &
                                       qsv_group(ibase, jbase, kbase) /= id_qsv_group_init) then
                                  if (qsv_group(ibase, jbase, kbase) /= qsv_group(i, j, k)) then
                                    call s_merge_groups(qsv_group, qsv_group(ibase, jbase, kbase), qsv_group(i, j, k))
                                  end if

                              ! if (ibase, jbase, kbase) is not grouped and (i, j, k) is grouped
                              else if (qsv_group(i, j, k) == id_qsv_group_init .and.  &
                                       qsv_group(ibase, jbase, kbase) /= id_qsv_group_init) then
                                  qsv_group(i, j, k) = qsv_group(ibase, jbase, kbase)

                              ! if (ibase, jbase, kbase) is grouped and (i, j, k) is not grouped
                              else
                                  qsv_group(ibase, jbase, kbase) = qsv_group(i, j, k)
                              end if
                          end if
                        end if
                      end do
                    end do
                  end do
                end if
              end do
            end do
          end if
        end do

        ! Merge groups across processors
        if (proc_rank == 0) print *, "Merge groups across processors"
        qsv_merge_x = .false.
        qsv_merge_z = .false.
        do l = 0, num_procs - 1
          ! Add paddings
          call s_add_paddings_integer(qsv_group, 1, qsv_group_padded)
          ! Get adjecent proc_rank
          call s_get_adjacent_proc_rank(1, pr_neg, pr_pos)
          if (proc_rank == l) then
            ! Merge qsv groups at the boundary in x-direction
            do k = 0, p
              do j = 0, n
                breakpoint = .false.
                do kk = k - 1, k + 1
                  do jj = j - 1, j + 1
                    if (qsv_group_padded(-1, jj, kk) /= pr_neg .and. &
                        qsv_group(0, j, k) /= proc_rank .and. &
                        qsv_group_padded(-1, jj, kk) /= qsv_group(0, j, k)) then
                      if (proc_rank_x == 0) then
                        where (qsv_group == qsv_group(0, j, k)) qsv_merge_x = .true.
                      end if
                      where (qsv_group == qsv_group(0, j, k)) qsv_group = qsv_group_padded(-1, jj, kk)
                      breakpoint = .true.
                    end if
                    if (breakpoint) exit
                  end do
                  if (breakpoint) exit
                end do
              end do
            end do
          end if

          ! Add paddings
          call s_add_paddings_integer(qsv_group, 1, qsv_group_padded)
          call s_add_paddings_logical(qsv_merge_x, 1, qsv_merge_x_padded)
          ! Get adjecent proc_rank
          call s_get_adjacent_proc_rank(2, pr_neg, pr_pos)
          if (proc_rank == l) then
            ! Merge qsv groups at the boundary in y-direction
            do k = 0, p
              do i = 0, m
                breakpoint = .false.
                do kk = k - 1, k + 1
                  do ii = i - 1, i + 1
                    if (qsv_group_padded(ii, -1, kk) /= pr_neg .and. &
                        qsv_group(i, 0, k) /= proc_rank .and. &
                        qsv_group_padded(ii, -1, kk) /= qsv_group(i, 0, k)) then
                      if (qsv_merge_x_padded(ii, -1, kk)) then
                        where (qsv_group == qsv_group(i, 0, k)) qsv_merge_x = .true.
                      end if
                      ! Have to deal with situation where two groups are separate in proc1 but they are connected by a point in proc2
                      WHERE (qsv_group == qsv_group(i, 0, k)) qsv_group = qsv_group_padded(ii, -1, kk)
                      breakpoint = .true.
                    end if
                    if (breakpoint) exit
                  end do
                  if (breakpoint) exit
                end do
              end do
            end do
          end if

          ! Add paddings
          call s_add_paddings_integer(qsv_group, 1, qsv_group_padded)
          call s_add_paddings_logical(qsv_merge_x, 1, qsv_merge_x_padded)
          ! Get adjecent proc_rank
          call s_get_adjacent_proc_rank(3, pr_neg, pr_pos)
          if (proc_rank == l) then
            ! Merge qsv groups at the boundary in z-direction
            do j = 0, n
              do i = 0, m
                breakpoint = .false.
                do jj = j - 1, j + 1
                  do ii = i - 1, i + 1
                    if (qsv_group_padded(ii, jj, -1) /= pr_neg .and. &
                        qsv_group(i, j, 0) /= proc_rank .and. &
                        qsv_group_padded(ii, jj, -1) /= qsv_group(i, j, 0)) then
                      if (qsv_merge_x_padded(ii, jj, -1)) then
                        where (qsv_group == qsv_group(i, 0, k)) qsv_merge_x = .true.
                      end if
                      if (proc_rank_z == 0) then
                        where (qsv_group == qsv_group(i, j, 0)) qsv_merge_z = .true.
                      end if
                      where (qsv_group == qsv_group(i, j, 0)) qsv_group = qsv_group_padded(ii, jj, -1)
                      breakpoint = .true.
                    end if
                    if (breakpoint) exit
                  end do
                  if (breakpoint) exit
                end do
              end do
            end do
          end if
        end do

        ! Set non-grouped points to group 0
        if (proc_rank == 0) print *, "Set non-grouped points to group 0"
        where (qsv_group == id_qsv_group_init) qsv_group = 0

    end subroutine s_identify_connected_groups

    impure subroutine s_merge_groups(qsv_group, id_qsv_group1, id_qsv_group2)
        integer, intent(in) :: id_qsv_group1, id_qsv_group2
        integer, dimension(0:m, 0:n, 0:p), intent(inout) :: qsv_group
        integer :: qsv_group_merged

        qsv_group_merged = min(id_qsv_group1, id_qsv_group2)
        where (qsv_group == id_qsv_group1) qsv_group = qsv_group_merged
        where (qsv_group == id_qsv_group2) qsv_group = qsv_group_merged

    end subroutine s_merge_groups

    impure subroutine s_remove_isolated_single_point(ibase, jbase, kbase, qsv_flag_in, qsv_flag_out)
        integer, intent(in) :: ibase, jbase, kbase
        logical, dimension(0:m, 0:n, 0:p), intent(in) :: qsv_flag_in
        logical, dimension(0:m, 0:n, 0:p), intent(out) :: qsv_flag_out
        integer :: i, j, k

        do k = kbase - 1, kbase + 1
          do j = jbase - 1, jbase + 1
            do i = ibase - 1, ibase + 1
              if (i >= 0 .and. i <= m .and. &
                  j >= 0 .and. j <= n .and. &
                  k >= 0 .and. k <= p) then
                if (.not. (i == ibase .and. j == jbase .and. k == kbase) &
                    .and. qsv_flag_in(i, j, k)) then
                  qsv_flag_out(ibase, jbase, kbase) = .true.
                  return
                end if
              end if
            end do
          end do
        end do
    end subroutine s_remove_isolated_single_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    impure subroutine s_apply_gaussian_filter(field_in, field_filtered, filter_size, y_idx_beg, y_idx_end)        
        real(wp), dimension(0:m, 0:n, 0:p), intent(in) :: field_in
        real(wp), dimension(0:m, 0:n, 0:p), intent(out) :: field_filtered
        real(wp), intent(in) :: filter_size
        integer, intent(in) :: y_idx_beg, y_idx_end

        integer :: pad_size
        real(wp), dimension(:, :, :), allocatable :: field
        real(wp), dimension(:), allocatable :: kernel_y, kernel_z

        real(wp) :: sigma, norm, dist_y, dist_z
        integer :: i, j, k, l, q, r, il, jq, kr, ii, jj, kk

        ! Compute pad size
        pad_size = nint((y_idx_end - y_idx_beg)/2._wp)
        if (pad_size > max(m + 1, n + 1)) then 
            print *, "pad_size is too large: ", pad_size, m + 1, n + 1
            call s_mpi_abort()
        end if

        ! Kernel parameters
        sigma = filter_size/20._wp

        ! Allocate memory
        allocate (field(0:m, -pad_size:n + pad_size, -pad_size:p + pad_size))
        allocate (kernel_y(-pad_size:pad_size))
        allocate (kernel_z(-pad_size:pad_size))
        
        ! Add paddings to field for each processor
        call s_add_paddings_real(field_in, 0, pad_size, pad_size, field)

        ! Filtering loops
        field_filtered = 0._wp
        do j = 0, n
            if (y_cc(j) >= y_cc_glb(y_idx_beg) .and. y_cc(j) <= y_cc_glb(y_idx_end)) then
                jj = j + proc_rank_y*(n + 1)
                do i = 0, m
                    do k = 0, p
                        kk = k + proc_rank_z*(p + 1)
                        ! Compute kernel
                        kernel_y = 0._wp
                        kernel_z = 0._wp
                        do l = -pad_size, pad_size
                            ! periodic in z
                            if (kk + l < 1) then
                                dist_z = (z_cc_glb(kk + l + p_glb) - (z_cb_glb(p_glb) - z_cb_glb(-1))) - z_cc_glb(kk)
                            else if (kk + l > p_glb) then
                                dist_z = (z_cc_glb(kk + l - p_glb) + (z_cb_glb(p_glb) - z_cb_glb(-1))) - z_cc_glb(kk)
                            else
                                dist_z = z_cc_glb(kk + l) - z_cc_glb(kk)
                            end if
                            ! kernel function
                            kernel_z(l) = exp(-0.5_wp * (dist_z / sigma)**2._wp)

                            ! non-periodic in y
                            if (jj + l < 1 .or. jj + l > n_glb) then
                                kernel_y(l) = 0._wp
                            else
                                dist_y = y_cc_glb(jj + l) - y_cc_glb(jj)
                                kernel_y(l) = exp(-0.5_wp * (dist_y / sigma)**2._wp)
                            end if
                        end do
                        kernel_y = kernel_y / sum(kernel_y)
                        kernel_z = kernel_z / sum(kernel_z)

                        ! Compute filtered field
                        do r = -pad_size, pad_size
                            do q = -pad_size, pad_size
                                jq = j + q
                                kr = k + r
                                field_filtered(i, j, k) = field_filtered(i, j, k) + field(i, jq, kr) * (kernel_y(q) * kernel_z(r))
                            end do
                        end do
                    end do
                end do
            end if
        end do

        deallocate (field, kernel_y, kernel_z)
    end subroutine s_apply_gaussian_filter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    impure subroutine s_add_paddings_xy(field_in, pad_size, field_out)
        real(wp), dimension(0:m, 0:n, 0:p), intent(in) :: field_in
        integer, intent(in) :: pad_size
        real(wp), dimension(-pad_size:m + pad_size, &
                            -pad_size:n + pad_size, &
                                    0:p), intent(out) :: field_out
        integer :: type_slice, pr_pos, pr_neg
        integer :: start_send, start_recv
        integer :: ierr

        ! Initialize output array
        field_out(0:m, 0:n, 0:p) = field_in(0:m, 0:n, 0:p)

        ! x-direction
        if (num_procs_x > 1) then
            ! Compute left neighbor
            if (proc_rank_x == 0) then
                call s_get_proc_rank(num_procs_x - 1, proc_rank_y, proc_rank_z, pr_neg)
            else
                call s_get_proc_rank(proc_rank_x - 1, proc_rank_y, proc_rank_z, pr_neg)
            end if
            ! Compute right neighbor
            if (proc_rank_x == num_procs_x - 1) then
                call s_get_proc_rank(0, proc_rank_y, proc_rank_z, pr_pos)
            else
                call s_get_proc_rank(proc_rank_x + 1, proc_rank_y, proc_rank_z, pr_pos)
            end if

            ! Create MPI_TYPE
            call create_slice_x_type(m + 1 + 2*pad_size, n + 1 + 2*pad_size, p + 1, pad_size, type_slice)

            ! Left padding
            start_send = m + 1 - pad_size
            start_recv = -pad_size
            call MPI_SENDRECV(field_out(start_send, -pad_size, 0), 1, type_slice, pr_pos, 0, &
                              field_out(start_recv, -pad_size, 0), 1, type_slice, pr_neg, 0, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

            ! Right padding
            start_send = 0
            start_recv = m + 1
            call MPI_SENDRECV(field_out(start_send, -pad_size, 0), 1, type_slice, pr_neg, 0, &
                              field_out(start_recv, -pad_size, 0), 1, type_slice, pr_pos, 0, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            
            ! Free MPI_TYPE
            call MPI_TYPE_FREE(type_slice, ierr)
        else
            ! Left padding
            field_out(-pad_size:-1, 0:n, 0:p) = field_in(m - pad_size + 1:m, 0:n, 0:p)
            ! Right padding
            field_out(m + 1:m + pad_size, 0:n, 0:p) = field_in(0:pad_size - 1, 0:n, 0:p)
        end if

        ! y-direction
        if (num_procs_y > 1) then
            ! Compute bottom neighbor
            if (proc_rank_y == 0) then
                call s_get_proc_rank(proc_rank_x, num_procs_y - 1, proc_rank_z, pr_neg)
            else
                call s_get_proc_rank(proc_rank_x, proc_rank_y - 1, proc_rank_z, pr_neg)
            end if
            ! Compute top neighbor
            if (proc_rank_y == num_procs_y - 1) then
                call s_get_proc_rank(proc_rank_x, 0, proc_rank_z, pr_pos)
            else
                call s_get_proc_rank(proc_rank_x, proc_rank_y + 1, proc_rank_z, pr_pos)
            end if

            ! Create MPI_TYPE
            call create_slice_y_type(m + 1 + 2*pad_size, n + 1 + 2*pad_size, p + 1, pad_size, type_slice)

            ! Bottom padding
            start_send = n + 1 - pad_size
            start_recv = -pad_size
            call MPI_SENDRECV(field_out(-pad_size, start_send, 0), 1, type_slice, pr_pos, 0, &
                              field_out(-pad_size, start_recv, 0), 1, type_slice, pr_neg, 0, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

            ! Top padding
            start_send = 0
            start_recv = n + 1
            call MPI_SENDRECV(field_out(-pad_size, start_send, 0), 1, type_slice, pr_neg, 0, &
                              field_out(-pad_size, start_recv, 0), 1, type_slice, pr_pos, 0, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

            ! Free MPI_TYPE
            call MPI_TYPE_FREE(type_slice, ierr)
        else
            ! Bottom padding
            field_out(0:m, -pad_size:-1, 0:p) = field_in(0:m, n - pad_size + 1:n, 0:p)
            ! Top padding
            field_out(0:m, n + 1:n + pad_size, 0:p) = field_in(0:m, 0:pad_size - 1, 0:p)
        end if

    contains
        subroutine create_slice_x_type(nx, ny, nz, len_slice, type_slice)
            integer, intent(in) :: nx, ny, nz, len_slice
            integer, intent(out) :: type_slice
            integer :: ierr
            call MPI_TYPE_VECTOR(ny*nz, len_slice, nx, mpi_p, type_slice, ierr)
            call MPI_TYPE_COMMIT(type_slice, ierr)
        end subroutine create_slice_x_type

        subroutine create_slice_y_type(nx, ny, nz, len_slice, type_slice)
            integer, intent(in) :: nx, ny, nz, len_slice
            integer, intent(out) :: type_slice
            integer :: ierr
            call MPI_TYPE_VECTOR(nz, len_slice*nx, nx*ny, mpi_p, type_slice, ierr)
            call MPI_TYPE_COMMIT(type_slice, ierr)
        end subroutine create_slice_y_type
    end subroutine s_add_paddings_xy

    impure subroutine s_add_paddings_real(field_in, pad_size_x, pad_size_y, pad_size_z, field_out)
        real(wp), dimension(0:m, 0:n, 0:p), intent(in) :: field_in
        integer, intent(in) :: pad_size_x, pad_size_y, pad_size_z
        real(wp), dimension(-pad_size_x:m + pad_size_x, &
                            -pad_size_y:n + pad_size_y, &
                            -pad_size_z:p + pad_size_z), intent(out) :: field_out
        integer :: type_slice, pr_pos, pr_neg
        integer :: start_send, start_recv
        integer :: ierr

        ! Initialize output array
        field_out(0:m, 0:n, 0:p) = field_in(0:m, 0:n, 0:p)

        ! x-direction
        if (pad_size_x > 0) then
            if (num_procs_x > 1) then
                ! Compute left neighbor
                if (proc_rank_x == 0) then
                    call s_get_proc_rank(num_procs_x - 1, proc_rank_y, proc_rank_z, pr_neg)
                else
                    call s_get_proc_rank(proc_rank_x - 1, proc_rank_y, proc_rank_z, pr_neg)
                end if
                ! Compute right neighbor
                if (proc_rank_x == num_procs_x - 1) then
                    call s_get_proc_rank(0, proc_rank_y, proc_rank_z, pr_pos)
                else
                    call s_get_proc_rank(proc_rank_x + 1, proc_rank_y, proc_rank_z, pr_pos)
                end if

                ! Create MPI_TYPE
                call create_slice_x_type(m + 1 + 2*pad_size_x, n + 1 + 2*pad_size_y, p + 1 + 2*pad_size_z, pad_size_x, type_slice)

                ! Left padding
                start_send = m + 1 - pad_size_x
                start_recv = -pad_size_x
                call MPI_SENDRECV(field_out(start_send, -pad_size_y, -pad_size_z), 1, type_slice, pr_pos, 0, &
                                field_out(start_recv, -pad_size_y, -pad_size_z), 1, type_slice, pr_neg, 0, &
                                MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                ! Right padding
                start_send = 0
                start_recv = m + 1
                call MPI_SENDRECV(field_out(start_send, -pad_size_y, -pad_size_z), 1, type_slice, pr_neg, 0, &
                                field_out(start_recv, -pad_size_y, -pad_size_z), 1, type_slice, pr_pos, 0, &
                                MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                
                ! Free MPI_TYPE
                call MPI_TYPE_FREE(type_slice, ierr)
            else
                ! Left padding
                field_out(-pad_size_x:-1, 0:n, 0:p) = field_in(m - pad_size_x + 1:m, 0:n, 0:p)
                ! Right padding
                field_out(m + 1:m + pad_size_x, 0:n, 0:p) = field_in(0:pad_size_x - 1, 0:n, 0:p)
            end if
        end if

        ! y-direction
        if (pad_size_y > 0) then
            if (num_procs_y > 1) then
                ! Compute bottom neighbor
                if (proc_rank_y == 0) then
                    call s_get_proc_rank(proc_rank_x, num_procs_y - 1, proc_rank_z, pr_neg)
                else
                    call s_get_proc_rank(proc_rank_x, proc_rank_y - 1, proc_rank_z, pr_neg)
                end if
                ! Compute top neighbor
                if (proc_rank_y == num_procs_y - 1) then
                    call s_get_proc_rank(proc_rank_x, 0, proc_rank_z, pr_pos)
                else
                    call s_get_proc_rank(proc_rank_x, proc_rank_y + 1, proc_rank_z, pr_pos)
                end if

                ! Create MPI_TYPE
                call create_slice_y_type(m + 1 + 2*pad_size_x, n + 1 + 2*pad_size_y, p + 1 + 2*pad_size_z, pad_size_y, type_slice)

                ! Bottom padding
                start_send = n + 1 - pad_size_y
                start_recv = -pad_size_y
                call MPI_SENDRECV(field_out(-pad_size_x, start_send, -pad_size_z), 1, type_slice, pr_pos, 0, &
                                field_out(-pad_size_x, start_recv, -pad_size_z), 1, type_slice, pr_neg, 0, &
                                MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                ! Top padding
                start_send = 0
                start_recv = n + 1
                call MPI_SENDRECV(field_out(-pad_size_x, start_send, -pad_size_z), 1, type_slice, pr_neg, 0, &
                                field_out(-pad_size_x, start_recv, -pad_size_z), 1, type_slice, pr_pos, 0, &
                                MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                ! Free MPI_TYPE
                call MPI_TYPE_FREE(type_slice, ierr)
            else
                ! Bottom padding
                field_out(0:m, -pad_size_y:-1, 0:p) = field_in(0:m, n - pad_size_y + 1:n, 0:p)
                ! Top padding
                field_out(0:m, n + 1:n + pad_size_y, 0:p) = field_in(0:m, 0:pad_size_y - 1, 0:p)
            end if
        end if

        ! z-direction
        if (pad_size_z > 0) then
            if (num_procs_z > 1) then
                ! Compute front neighbor
                if (proc_rank_z == 0) then
                    call s_get_proc_rank(proc_rank_x, proc_rank_y, num_procs_z - 1, pr_neg)
                else
                    call s_get_proc_rank(proc_rank_x, proc_rank_y, proc_rank_z - 1, pr_neg)
                end if
                ! Compute back neighbor
                if (proc_rank_z == num_procs_z - 1) then
                    call s_get_proc_rank(proc_rank_x, proc_rank_y, 0, pr_pos)
                else
                    call s_get_proc_rank(proc_rank_x, proc_rank_y, proc_rank_z + 1, pr_pos)
                end if

                ! Create MPI_TYPE
                call create_slice_z_type(m + 1 + 2*pad_size_x, n + 1 + 2*pad_size_y, p + 1 + 2*pad_size_z, pad_size_z, type_slice)

                ! Bottom padding
                start_send = p + 1 - pad_size_z
                start_recv = -pad_size_z
                call MPI_SENDRECV(field_out(-pad_size_x, -pad_size_y, start_send), 1, type_slice, pr_pos, 0, &
                                field_out(-pad_size_x, -pad_size_y, start_recv), 1, type_slice, pr_neg, 0, &
                                MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                ! Top padding
                start_send = 0
                start_recv = p + 1
                call MPI_SENDRECV(field_out(-pad_size_x, -pad_size_y, start_send), 1, type_slice, pr_neg, 0, &
                                field_out(-pad_size_x, -pad_size_y, start_recv), 1, type_slice, pr_pos, 0, &
                                MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                ! Free MPI_TYPE
                call MPI_TYPE_FREE(type_slice, ierr)
            else
                ! Front padding
                field_out(0:m, 0:p, -pad_size_z:-1) = field_in(0:m, 0:n, p - pad_size_z + 1:p)
                ! Back padding
                field_out(0:m, 0:n, p + 1:p + pad_size_z) = field_in(0:m, 0:n, 0:pad_size_z - 1)
            end if
        end if

    contains
        subroutine create_slice_x_type(nx, ny, nz, len_slice, type_slice)
            integer, intent(in) :: nx, ny, nz, len_slice
            integer, intent(out) :: type_slice
            integer :: ierr
            call MPI_TYPE_VECTOR(ny*nz, len_slice, nx, mpi_p, type_slice, ierr)
            call MPI_TYPE_COMMIT(type_slice, ierr)
        end subroutine create_slice_x_type

        subroutine create_slice_y_type(nx, ny, nz, len_slice, type_slice)
            integer, intent(in) :: nx, ny, nz, len_slice
            integer, intent(out) :: type_slice
            integer :: ierr
            call MPI_TYPE_VECTOR(nz, len_slice*nx, nx*ny, mpi_p, type_slice, ierr)
            call MPI_TYPE_COMMIT(type_slice, ierr)
        end subroutine create_slice_y_type

        subroutine create_slice_z_type(nx, ny, nz, len_slice, type_slice)
            integer, intent(in) :: nx, ny, nz, len_slice
            integer, intent(out) :: type_slice
            integer :: ierr
            call MPI_TYPE_CONTIGUOUS(nx*ny*len_slice, mpi_p, type_slice, ierr)
            call MPI_TYPE_COMMIT(type_slice, ierr)
        end subroutine create_slice_z_type
    end subroutine s_add_paddings_real

    impure subroutine s_add_paddings_integer(field_in, pad_size, field_out)
        integer, dimension(0:m, 0:n, 0:p), intent(in) :: field_in
        integer, intent(in) :: pad_size
        integer,  dimension(-pad_size:m + pad_size, &
                            -pad_size:n + pad_size, &
                            -pad_size:p + pad_size), intent(out) :: field_out
        integer :: type_slice, pr_pos, pr_neg
        integer :: start_send, start_recv
        integer :: ierr

        ! Initialize output array
        field_out(0:m, 0:n, 0:p) = field_in(0:m, 0:n, 0:p)

        ! x-direction
        if (num_procs_x > 1) then
            ! Compute left neighbor
            if (proc_rank_x == 0) then
                call s_get_proc_rank(num_procs_x - 1, proc_rank_y, proc_rank_z, pr_neg)
            else
                call s_get_proc_rank(proc_rank_x - 1, proc_rank_y, proc_rank_z, pr_neg)
            end if
            ! Compute right neighbor
            if (proc_rank_x == num_procs_x - 1) then
                call s_get_proc_rank(0, proc_rank_y, proc_rank_z, pr_pos)
            else
                call s_get_proc_rank(proc_rank_x + 1, proc_rank_y, proc_rank_z, pr_pos)
            end if

            ! Create MPI_TYPE
            call create_slice_x_type(m + 1 + 2*pad_size, n + 1 + 2*pad_size, p + 1 + 2*pad_size, pad_size, type_slice)

            ! Left padding
            start_send = m + 1 - pad_size
            start_recv = -pad_size
            call MPI_SENDRECV(field_out(start_send, -pad_size, -pad_size), 1, type_slice, pr_pos, 0, &
                              field_out(start_recv, -pad_size, -pad_size), 1, type_slice, pr_neg, 0, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

            ! Right padding
            start_send = 0
            start_recv = m + 1
            call MPI_SENDRECV(field_out(start_send, -pad_size, -pad_size), 1, type_slice, pr_neg, 0, &
                              field_out(start_recv, -pad_size, -pad_size), 1, type_slice, pr_pos, 0, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            
            ! Free MPI_TYPE
            call MPI_TYPE_FREE(type_slice, ierr)
        else
            ! Left padding
            field_out(-pad_size:-1, 0:n, 0:p) = field_in(m - pad_size + 1:m, 0:n, 0:p)
            ! Right padding
            field_out(m + 1:m + pad_size, 0:n, 0:p) = field_in(0:pad_size - 1, 0:n, 0:p)
        end if

        ! y-direction
        if (num_procs_y > 1) then
            ! Compute bottom neighbor
            if (proc_rank_y == 0) then
                call s_get_proc_rank(proc_rank_x, num_procs_y - 1, proc_rank_z, pr_neg)
            else
                call s_get_proc_rank(proc_rank_x, proc_rank_y - 1, proc_rank_z, pr_neg)
            end if
            ! Compute top neighbor
            if (proc_rank_y == num_procs_y - 1) then
                call s_get_proc_rank(proc_rank_x, 0, proc_rank_z, pr_pos)
            else
                call s_get_proc_rank(proc_rank_x, proc_rank_y + 1, proc_rank_z, pr_pos)
            end if

            ! Create MPI_TYPE
            call create_slice_y_type(m + 1 + 2*pad_size, n + 1 + 2*pad_size, p + 1 + 2*pad_size, pad_size, type_slice)

            ! Bottom padding
            start_send = n + 1 - pad_size
            start_recv = -pad_size
            call MPI_SENDRECV(field_out(-pad_size, start_send, -pad_size), 1, type_slice, pr_pos, 0, &
                              field_out(-pad_size, start_recv, -pad_size), 1, type_slice, pr_neg, 0, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

            ! Top padding
            start_send = 0
            start_recv = n + 1
            call MPI_SENDRECV(field_out(-pad_size, start_send, -pad_size), 1, type_slice, pr_neg, 0, &
                              field_out(-pad_size, start_recv, -pad_size), 1, type_slice, pr_pos, 0, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

            ! Free MPI_TYPE
            call MPI_TYPE_FREE(type_slice, ierr)
        else
            ! Bottom padding
            field_out(0:m, -pad_size:-1, 0:p) = field_in(0:m, n - pad_size + 1:n, 0:p)
            ! Top padding
            field_out(0:m, n + 1:n + pad_size, 0:p) = field_in(0:m, 0:pad_size - 1, 0:p)
        end if

        ! z-direction
        if (num_procs_z > 1) then
            ! Compute front neighbor
            if (proc_rank_z == 0) then
                call s_get_proc_rank(proc_rank_x, proc_rank_y, num_procs_z - 1, pr_neg)
            else
                call s_get_proc_rank(proc_rank_x, proc_rank_y, proc_rank_z - 1, pr_neg)
            end if
            ! Compute back neighbor
            if (proc_rank_z == num_procs_z - 1) then
                call s_get_proc_rank(proc_rank_x, proc_rank_y, 0, pr_pos)
            else
                call s_get_proc_rank(proc_rank_x, proc_rank_y, proc_rank_z + 1, pr_pos)
            end if

            ! Create MPI_TYPE
            call create_slice_z_type(m + 1 + 2*pad_size, n + 1 + 2*pad_size, p + 1 + 2*pad_size, pad_size, type_slice)

            ! Bottom padding
            start_send = p + 1 - pad_size
            start_recv = -pad_size
            call MPI_SENDRECV(field_out(-pad_size, -pad_size, start_send), 1, type_slice, pr_pos, 0, &
                              field_out(-pad_size, -pad_size, start_recv), 1, type_slice, pr_neg, 0, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

            ! Top padding
            start_send = 0
            start_recv = p + 1
            call MPI_SENDRECV(field_out(-pad_size, -pad_size, start_send), 1, type_slice, pr_neg, 0, &
                              field_out(-pad_size, -pad_size, start_recv), 1, type_slice, pr_pos, 0, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

            ! Free MPI_TYPE
            call MPI_TYPE_FREE(type_slice, ierr)
        else
            ! Front padding
            field_out(0:m, 0:p, -pad_size:-1) = field_in(0:m, 0:n, p - pad_size + 1:p)
            ! Back padding
            field_out(0:m, 0:n, p + 1:p + pad_size) = field_in(0:m, 0:n, 0:pad_size - 1)
        end if

    contains
        subroutine create_slice_x_type(nx, ny, nz, len_slice, type_slice)
            integer, intent(in) :: nx, ny, nz, len_slice
            integer, intent(out) :: type_slice
            integer :: ierr
            call MPI_TYPE_VECTOR(ny*nz, len_slice, nx, mpi_integer, type_slice, ierr)
            call MPI_TYPE_COMMIT(type_slice, ierr)
        end subroutine create_slice_x_type

        subroutine create_slice_y_type(nx, ny, nz, len_slice, type_slice)
            integer, intent(in) :: nx, ny, nz, len_slice
            integer, intent(out) :: type_slice
            integer :: ierr
            call MPI_TYPE_VECTOR(nz, len_slice*nx, nx*ny, mpi_integer, type_slice, ierr)
            call MPI_TYPE_COMMIT(type_slice, ierr)
        end subroutine create_slice_y_type

        subroutine create_slice_z_type(nx, ny, nz, len_slice, type_slice)
            integer, intent(in) :: nx, ny, nz, len_slice
            integer, intent(out) :: type_slice
            integer :: ierr
            call MPI_TYPE_CONTIGUOUS(nx*ny*len_slice, mpi_integer, type_slice, ierr)
            call MPI_TYPE_COMMIT(type_slice, ierr)
        end subroutine create_slice_z_type
    end subroutine s_add_paddings_integer

    impure subroutine s_add_paddings_logical(field_in, pad_size, field_out)
        logical, dimension(0:m, 0:n, 0:p), intent(in) :: field_in
        integer, intent(in) :: pad_size
        logical,  dimension(-pad_size:m + pad_size, &
                            -pad_size:n + pad_size, &
                            -pad_size:p + pad_size), intent(out) :: field_out
        integer :: type_slice, pr_pos, pr_neg
        integer :: start_send, start_recv
        integer :: ierr

        ! Initialize output array
        field_out(0:m, 0:n, 0:p) = field_in(0:m, 0:n, 0:p)

        ! x-direction
        if (num_procs_x > 1) then
            ! Compute left neighbor
            if (proc_rank_x == 0) then
                call s_get_proc_rank(num_procs_x - 1, proc_rank_y, proc_rank_z, pr_neg)
            else
                call s_get_proc_rank(proc_rank_x - 1, proc_rank_y, proc_rank_z, pr_neg)
            end if
            ! Compute right neighbor
            if (proc_rank_x == num_procs_x - 1) then
                call s_get_proc_rank(0, proc_rank_y, proc_rank_z, pr_pos)
            else
                call s_get_proc_rank(proc_rank_x + 1, proc_rank_y, proc_rank_z, pr_pos)
            end if

            ! Create MPI_TYPE
            call create_slice_x_type(m + 1 + 2*pad_size, n + 1 + 2*pad_size, p + 1 + 2*pad_size, pad_size, type_slice)

            ! Left padding
            start_send = m + 1 - pad_size
            start_recv = -pad_size
            call MPI_SENDRECV(field_out(start_send, -pad_size, -pad_size), 1, type_slice, pr_pos, 0, &
                              field_out(start_recv, -pad_size, -pad_size), 1, type_slice, pr_neg, 0, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

            ! Right padding
            start_send = 0
            start_recv = m + 1
            call MPI_SENDRECV(field_out(start_send, -pad_size, -pad_size), 1, type_slice, pr_neg, 0, &
                              field_out(start_recv, -pad_size, -pad_size), 1, type_slice, pr_pos, 0, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            
            ! Free MPI_TYPE
            call MPI_TYPE_FREE(type_slice, ierr)
        else
            ! Left padding
            field_out(-pad_size:-1, 0:n, 0:p) = field_in(m - pad_size + 1:m, 0:n, 0:p)
            ! Right padding
            field_out(m + 1:m + pad_size, 0:n, 0:p) = field_in(0:pad_size - 1, 0:n, 0:p)
        end if

        ! y-direction
        if (num_procs_y > 1) then
            ! Compute bottom neighbor
            if (proc_rank_y == 0) then
                call s_get_proc_rank(proc_rank_x, num_procs_y - 1, proc_rank_z, pr_neg)
            else
                call s_get_proc_rank(proc_rank_x, proc_rank_y - 1, proc_rank_z, pr_neg)
            end if
            ! Compute top neighbor
            if (proc_rank_y == num_procs_y - 1) then
                call s_get_proc_rank(proc_rank_x, 0, proc_rank_z, pr_pos)
            else
                call s_get_proc_rank(proc_rank_x, proc_rank_y + 1, proc_rank_z, pr_pos)
            end if

            ! Create MPI_TYPE
            call create_slice_y_type(m + 1 + 2*pad_size, n + 1 + 2*pad_size, p + 1 + 2*pad_size, pad_size, type_slice)

            ! Bottom padding
            start_send = n + 1 - pad_size
            start_recv = -pad_size
            call MPI_SENDRECV(field_out(-pad_size, start_send, -pad_size), 1, type_slice, pr_pos, 0, &
                              field_out(-pad_size, start_recv, -pad_size), 1, type_slice, pr_neg, 0, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

            ! Top padding
            start_send = 0
            start_recv = n + 1
            call MPI_SENDRECV(field_out(-pad_size, start_send, -pad_size), 1, type_slice, pr_neg, 0, &
                              field_out(-pad_size, start_recv, -pad_size), 1, type_slice, pr_pos, 0, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

            ! Free MPI_TYPE
            call MPI_TYPE_FREE(type_slice, ierr)
        else
            ! Bottom padding
            field_out(0:m, -pad_size:-1, 0:p) = field_in(0:m, n - pad_size + 1:n, 0:p)
            ! Top padding
            field_out(0:m, n + 1:n + pad_size, 0:p) = field_in(0:m, 0:pad_size - 1, 0:p)
        end if

        ! z-direction
        if (num_procs_z > 1) then
            ! Compute front neighbor
            if (proc_rank_z == 0) then
                call s_get_proc_rank(proc_rank_x, proc_rank_y, num_procs_z - 1, pr_neg)
            else
                call s_get_proc_rank(proc_rank_x, proc_rank_y, proc_rank_z - 1, pr_neg)
            end if
            ! Compute back neighbor
            if (proc_rank_z == num_procs_z - 1) then
                call s_get_proc_rank(proc_rank_x, proc_rank_y, 0, pr_pos)
            else
                call s_get_proc_rank(proc_rank_x, proc_rank_y, proc_rank_z + 1, pr_pos)
            end if

            ! Create MPI_TYPE
            call create_slice_z_type(m + 1 + 2*pad_size, n + 1 + 2*pad_size, p + 1 + 2*pad_size, pad_size, type_slice)

            ! Bottom padding
            start_send = p + 1 - pad_size
            start_recv = -pad_size
            call MPI_SENDRECV(field_out(-pad_size, -pad_size, start_send), 1, type_slice, pr_pos, 0, &
                              field_out(-pad_size, -pad_size, start_recv), 1, type_slice, pr_neg, 0, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

            ! Top padding
            start_send = 0
            start_recv = p + 1
            call MPI_SENDRECV(field_out(-pad_size, -pad_size, start_send), 1, type_slice, pr_neg, 0, &
                              field_out(-pad_size, -pad_size, start_recv), 1, type_slice, pr_pos, 0, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

            ! Free MPI_TYPE
            call MPI_TYPE_FREE(type_slice, ierr)
        else
            ! Front padding
            field_out(0:m, 0:p, -pad_size:-1) = field_in(0:m, 0:n, p - pad_size + 1:p)
            ! Back padding
            field_out(0:m, 0:n, p + 1:p + pad_size) = field_in(0:m, 0:n, 0:pad_size - 1)
        end if

    contains
        subroutine create_slice_x_type(nx, ny, nz, len_slice, type_slice)
            integer, intent(in) :: nx, ny, nz, len_slice
            integer, intent(out) :: type_slice
            integer :: ierr
            call MPI_TYPE_VECTOR(ny*nz, len_slice, nx, mpi_logical, type_slice, ierr)
            call MPI_TYPE_COMMIT(type_slice, ierr)
        end subroutine create_slice_x_type

        subroutine create_slice_y_type(nx, ny, nz, len_slice, type_slice)
            integer, intent(in) :: nx, ny, nz, len_slice
            integer, intent(out) :: type_slice
            integer :: ierr
            call MPI_TYPE_VECTOR(nz, len_slice*nx, nx*ny, mpi_logical, type_slice, ierr)
            call MPI_TYPE_COMMIT(type_slice, ierr)
        end subroutine create_slice_y_type

        subroutine create_slice_z_type(nx, ny, nz, len_slice, type_slice)
            integer, intent(in) :: nx, ny, nz, len_slice
            integer, intent(out) :: type_slice
            integer :: ierr
            call MPI_TYPE_CONTIGUOUS(nx*ny*len_slice, mpi_logical, type_slice, ierr)
            call MPI_TYPE_COMMIT(type_slice, ierr)
        end subroutine create_slice_z_type
    end subroutine s_add_paddings_logical
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    impure subroutine s_get_proc_rank_xyz()
        integer :: proc_rank_tmp
        proc_rank_z = mod(proc_rank, num_procs_z)
        proc_rank_tmp = (proc_rank - proc_rank_z)/num_procs_z
        proc_rank_y = mod(proc_rank_tmp, num_procs_y)
        proc_rank_x = (proc_rank_tmp - proc_rank_y)/num_procs_y
    end subroutine s_get_proc_rank_xyz

    impure subroutine s_get_proc_rank(prx, pry, prz, pr)
        integer, intent(in) :: prx, pry, prz 
        integer, intent(out) :: pr

        pr = prx*num_procs_y*num_procs_z + pry*num_procs_z + prz
    end subroutine s_get_proc_rank

    impure subroutine s_get_adjacent_proc_rank(dir, pr_neg, pr_pos)
      integer, intent(in) :: dir
      integer, intent(out) :: pr_neg, pr_pos

      if (dir == 1) then
        if (num_procs_x > 1) then
          ! Compute left neighbor
          if (proc_rank_x == 0) then
              call s_get_proc_rank(num_procs_x - 1, proc_rank_y, proc_rank_z, pr_neg)
          else
              call s_get_proc_rank(proc_rank_x - 1, proc_rank_y, proc_rank_z, pr_neg)
          end if
          ! Compute right neighbor
          if (proc_rank_x == num_procs_x - 1) then
              call s_get_proc_rank(0, proc_rank_y, proc_rank_z, pr_pos)
          else
              call s_get_proc_rank(proc_rank_x + 1, proc_rank_y, proc_rank_z, pr_pos)
          end if
        else
          pr_neg = 0
          pr_pos = 0
        end if
      else if (dir == 2) then
        if (num_procs_y > 1) then
            ! Compute bottom neighbor
            if (proc_rank_y == 0) then
                call s_get_proc_rank(proc_rank_x, num_procs_y - 1, proc_rank_z, pr_neg)
            else
                call s_get_proc_rank(proc_rank_x, proc_rank_y - 1, proc_rank_z, pr_neg)
            end if
            ! Compute top neighbor
            if (proc_rank_y == num_procs_y - 1) then
                call s_get_proc_rank(proc_rank_x, 0, proc_rank_z, pr_pos)
            else
                call s_get_proc_rank(proc_rank_x, proc_rank_y + 1, proc_rank_z, pr_pos)
            end if
        else
          pr_neg = 0
          pr_pos = 0
        end if
      else if (dir == 3) then
        if (num_procs_z > 1) then
          ! Compute front neighbor
          if (proc_rank_z == 0) then
              call s_get_proc_rank(proc_rank_x, proc_rank_y, num_procs_z - 1, pr_neg)
          else
              call s_get_proc_rank(proc_rank_x, proc_rank_y, proc_rank_z - 1, pr_neg)
          end if
          ! Compute back neighbor
          if (proc_rank_z == num_procs_z - 1) then
              call s_get_proc_rank(proc_rank_x, proc_rank_y, 0, pr_pos)
          else
              call s_get_proc_rank(proc_rank_x, proc_rank_y, proc_rank_z + 1, pr_pos)
          end if
        else
          pr_neg = 0
          pr_pos = 0
        end if
      end if

    end subroutine s_get_adjacent_proc_rank
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !>  This subroutine gets as inputs the conservative variables
        !!      and density. From those inputs, it proceeds to calculate
        !!      the values of the numerical Schlieren function, which are
        !!      subsequently stored in the derived flow quantity storage
        !!      variable, q_sf.
        !!  @param q_cons_vf Conservative variables
        !!  @param q_sf Numerical Schlieren function
    impure subroutine s_derive_numerical_schlieren_function(q_cons_vf, q_sf)

        type(scalar_field), &
            dimension(sys_size), &
            intent(in) :: q_cons_vf

        real(wp), &
            dimension(-offset_x%beg:m + offset_x%end, &
                      -offset_y%beg:n + offset_y%end, &
                      -offset_z%beg:p + offset_z%end), &
            intent(inout) :: q_sf

        real(wp) :: drho_dx, drho_dy, drho_dz !<
            !! Spatial derivatives of the density in the x-, y- and z-directions

        real(wp), dimension(2) :: gm_rho_max !<
            !! Maximum value of the gradient magnitude (gm) of the density field
            !! in entire computational domain and not just the local sub-domain.
            !! The first position in the variable contains the maximum value and
            !! the second contains the rank of the processor on which it occurred.

        integer :: i, j, k, l !< Generic loop iterators

        ! Computing Gradient Magnitude of Density

        ! Contributions from the x- and y-coordinate directions
        do l = -offset_z%beg, p + offset_z%end
            do k = -offset_y%beg, n + offset_y%end
                do j = -offset_x%beg, m + offset_x%end

                    drho_dx = 0._wp
                    drho_dy = 0._wp

                    do i = -fd_number, fd_number
                        drho_dx = drho_dx + fd_coeff_x(i, j)*rho_sf(i + j, k, l)
                        drho_dy = drho_dy + fd_coeff_y(i, k)*rho_sf(j, i + k, l)
                    end do

                    gm_rho_sf(j, k, l) = drho_dx*drho_dx + drho_dy*drho_dy

                end do
            end do
        end do

        ! Contribution from the z-coordinate direction
        if (p > 0) then
            do l = -offset_z%beg, p + offset_z%end
                do k = -offset_y%beg, n + offset_y%end
                    do j = -offset_x%beg, m + offset_x%end

                        drho_dz = 0._wp

                        do i = -fd_number, fd_number
                            if (grid_geometry == 3) then
                                drho_dz = drho_dz + fd_coeff_z(i, l)/y_cc(k)* &
                                          rho_sf(j, k, i + l)
                            else
                                drho_dz = drho_dz + fd_coeff_z(i, l)* &
                                          rho_sf(j, k, i + l)
                            end if
                        end do

                        gm_rho_sf(j, k, l) = gm_rho_sf(j, k, l) &
                                             + drho_dz*drho_dz

                    end do
                end do
            end do
        end if

        ! Up until now, only the dot product of the gradient of the density
        ! field has been calculated and stored in the gradient magnitude of
        ! density variable. So now we proceed to take the square-root as to
        ! complete the desired calculation.
        gm_rho_sf = sqrt(gm_rho_sf)

        ! Determining the local maximum of the gradient magnitude of density
        ! and bookkeeping the result, along with rank of the local processor
        gm_rho_max = (/maxval(gm_rho_sf), real(proc_rank, wp)/)

        ! Comparing the local maximum gradient magnitude of the density on
        ! this processor to the those computed on the remaining processors.
        ! This allows for the global maximum to be computed and the rank of
        ! the processor on which it has occurred to be recorded.
        if (num_procs > 1) call s_mpi_reduce_maxloc(gm_rho_max)

        ! Computing Numerical Schlieren Function

        ! The form of the numerical Schlieren function depends on the choice
        ! of the multicomponent flow model. For the gamma/pi_inf model, the
        ! exponential of the negative, normalized, gradient magnitude of the
        ! density is computed. For the volume fraction model, the amplitude
        ! of the exponential's inside is also modulated with respect to the
        ! identity of the fluid in which the function is evaluated. For more
        ! information, refer to Marquina and Mulet (2003).

        if (model_eqns == 1) then                    ! Gamma/pi_inf model
            q_sf = -gm_rho_sf/gm_rho_max(1)

        else                                        ! Volume fraction model
            do l = -offset_z%beg, p + offset_z%end
                do k = -offset_y%beg, n + offset_y%end
                    do j = -offset_x%beg, m + offset_x%end

                        q_sf(j, k, l) = 0._wp

                        do i = 1, adv_idx%end - E_idx
                            q_sf(j, k, l) = &
                                q_sf(j, k, l) - schlieren_alpha(i)* &
                                q_cons_vf(i + E_idx)%sf(j, k, l)* &
                                gm_rho_sf(j, k, l)/gm_rho_max(1)
                        end do
                    end do
                end do
            end do
        end if

        ! Up until now, only the inside of the exponential of the numerical
        ! Schlieren function has been evaluated and stored. Then, to finish
        ! the computation, the exponential of the inside quantity is taken.
        q_sf = exp(q_sf)

    end subroutine s_derive_numerical_schlieren_function

    !>  Deallocation procedures for the module
    impure subroutine s_finalize_derived_variables_module

        ! Deallocating the variable containing the gradient magnitude of the
        ! density field provided that the numerical Schlieren function was
        ! was outputted during the post-process
        if (schlieren_wrt) deallocate (gm_rho_sf)

        ! Deallocating the variables that might have been used to bookkeep
        ! the finite-difference coefficients in the x-, y- and z-directions
        if (allocated(fd_coeff_x)) deallocate (fd_coeff_x)
        if (allocated(fd_coeff_y)) deallocate (fd_coeff_y)
        if (allocated(fd_coeff_z)) deallocate (fd_coeff_z)

    end subroutine s_finalize_derived_variables_module

end module m_derived_variables
