#!/usr/bin/env python3
import math
import json

## Fluid property
# water
gam_w = 6.12
pi_inf_w = 3.43e8
rho_w = 1000
mul0 = 1.002e-03  # viscosity
ss = 0.07275  # surface tension
pv = 2.3388e03  # vapor pressure

# vapor
gamma_v = 1.33
pi_inf_v = 0.0

# air
gam_a = 1.4
pi_inf_a = 0.0
rho_a_pre_shock = 1.17
rho_a_post_shock = 2.18

## Field
D = 0.048     # Droplet diameter
xbeg = -4*D   # left domain boundary
xend = 20*D   # right domain boundary
ybeg = 0.0    # bottom domain boundary
yend = 6*D    # top domain boundary
y_a = -5.7*D  # y-grid stretching parameter
y_b = 5.7*D   # y-grid stretching parameter
a_y = 3.67    # y-grid stretching parameter
loops_y = 2   # y-grid stretching parameter
R0ref = 10e-6 # Subgrid bubble radius reference scale

## Patch 1
p1_xc = 8*D
p1_yc = 6*D
p1_lx = 24*D
p1_ly = 14*D
p1_vel1 = 0.0
p1_vel2 = 0.0
p1_pres = 101325.0
p1_alpha1 = 1e-9
p1_alpha2 = 1 - p1_alpha1
p1_alpharho1 = p1_alpha1*rho_w
p1_alpharho2 = p1_alpha2*rho_a_pre_shock
p1_r0 = R0ref
p1_v0 = 0.0

## Patch 2
p2_xc = -2.5*D
p2_yc = 6*D
p2_lx = 3*D
p2_ly = 14*D
p2_vel1 = 226
p2_vel2 = 0.0
p2_pres = 238558.0
p2_alpha1 = 1e-9
p2_alpha2 = 1 - p2_alpha1
p2_alpharho1 = p2_alpha1*rho_w
p2_alpharho2 = p2_alpha2*rho_a_post_shock
p2_r0 = R0ref
p2_v0 = 0.0

## Patch 3
p3_xc = 0.0
p3_yc = 0.0
p3_radius = D/2
p3_vel1 = 0.0
p3_vel2 = 0.0
p3_pres = 101325.0
p3_alpha1 = 1 - 1e-9
p3_alpha2 = 1 - p3_alpha1
p3_alpharho1 = p3_alpha1*rho_w
p3_alpharho2 = p3_alpha2*rho_a_pre_shock
p3_r0 = R0ref
p3_v0 = 0.0

## REFERENCE VALUES
rho0 = rho_w
x0 = R0ref
p0 = 101325.0
u0 = math.sqrt(p0/rho0)

## Bubbles
icsg_vf = 1e-5
Ca = (p0 - pv)/p0
Web = rho0*u0**2*x0/ss
Re_inv = mul0/(rho0*u0*x0)

c_l = math.sqrt(1.4 * 238558.0 / 1)

## Grid
Ny = 299.0
Nx = 1199.0
dx = 0.25 / Nx  # 8.3e-6

time_end = 0.005  # 50us
cfl = 0.25

dt = cfl * dx / c_l  # 5.3E-9
Nt = int(time_end / dt)  # 10000

print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "F",
            # Computational Domain Parameters
            "x_domain%beg": xbeg/x0,
            "x_domain%end": xend/x0,
            "y_domain%beg": ybeg/x0,
            "y_domain%end": yend/x0,
            "stretch_y": "T",
            "a_y": a_y,
            "y_a": y_a/x0,
            "y_b": y_b/x0,
            "loops_y": loops_y,
            "m": int(Nx),
            "n": int(Ny),
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": Nt,
            # Simulation Algorithm Parameters
            "num_patches": 3,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            # "mpp_lim": "T",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -6,
            "bc_x%end": -6,
            "bc_y%beg": -2,
            "bc_y%end": -3,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1: Background
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": p1_xc/x0,
            "patch_icpp(1)%y_centroid": p1_yc/x0,
            "patch_icpp(1)%length_x": p1_lx/x0,
            "patch_icpp(1)%length_y": p1_ly/x0,
            "patch_icpp(1)%vel(1)": p1_vel1/u0,
            "patch_icpp(1)%vel(2)": p1_vel2/u0,
            "patch_icpp(1)%pres": p1_pres/p0,
            "patch_icpp(1)%alpha_rho(1)": p1_alpharho1/rho0,
            "patch_icpp(1)%alpha_rho(2)": p1_alpharho2/rho0,
            "patch_icpp(1)%alpha(1)": p1_alpha1,
            "patch_icpp(1)%alpha(2)": p1_alpha2,
            "patch_icpp(1)%r0": p1_r0/x0,
            "patch_icpp(1)%v0": p1_v0/u0,
            # Patch 2: Shocked state
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": p2_xc/x0,
            "patch_icpp(2)%y_centroid": p2_yc/x0,
            "patch_icpp(2)%length_x": p2_lx/x0,
            "patch_icpp(2)%length_y": p2_ly/x0,
            "patch_icpp(2)%vel(1)": p2_vel1/u0,
            "patch_icpp(2)%vel(2)": p2_vel2/u0,
            "patch_icpp(2)%pres": p2_pres/p0,
            "patch_icpp(2)%alpha_rho(1)": p2_alpharho1/rho0,
            "patch_icpp(2)%alpha_rho(2)": p2_alpharho2/rho0,
            "patch_icpp(2)%alpha(1)": p2_alpha1,
            "patch_icpp(2)%alpha(2)": p2_alpha2,
            "patch_icpp(2)%r0": p2_r0/x0,
            "patch_icpp(2)%v0": p2_v0/u0,
            # Patch 3: Bubble
            "patch_icpp(3)%geometry": 2,
            "patch_icpp(3)%x_centroid": p3_xc/x0,
            "patch_icpp(3)%y_centroid": p3_yc/x0,
            "patch_icpp(3)%radius": p3_radius/x0,
            "patch_icpp(3)%alter_patch(1)": "T",
            "patch_icpp(3)%vel(1)": p3_vel1/u0,
            "patch_icpp(3)%vel(2)": p3_vel2/u0,
            "patch_icpp(3)%pres": p3_pres/p0,
            "patch_icpp(3)%alpha_rho(1)": p3_alpharho1/rho0,
            "patch_icpp(3)%alpha_rho(2)": p3_alpharho2/rho0,
            "patch_icpp(3)%alpha(1)": p3_alpha1,
            "patch_icpp(3)%alpha(2)": p3_alpha2,
            "patch_icpp(3)%r0": p3_r0/x0,
            "patch_icpp(3)%v0": p3_v0/u0,
            # Subgrid bubble model
            "bubbles_euler": "T",
            "bubble_model": 2,
            "icsg": "T",
            "icsg_vf": icsg_vf,
            "icsg_patch": 3,
            "adv_n": "T",
            "polytropic": "T",
            "polydisperse": "F",
            "nb": 1,
            "Ca": Ca,
            "Web": Web,
            "Re_inv": Re_inv,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0 / (gam_w - 1.0),
            "fluid_pp(1)%pi_inf": gam_w * (pi_inf_w/p0) / (gam_w - 1.0),
            "fluid_pp(2)%gamma": 1.0 / (gam_a - 1.0),
            "fluid_pp(2)%pi_inf": gam_a * (pi_inf_a/p0) / (gam_a - 1.0),
            "fluid_pp(3)%gamma": 1.0 / (gam_a - 1.0),
            "fluid_pp(3)%pi_inf": gam_a * (pi_inf_a/p0) / (gam_a - 1.0),
        }
    )
)
