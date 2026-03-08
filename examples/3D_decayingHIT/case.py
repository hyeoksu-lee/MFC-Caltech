#!/usr/bin/env python3
import math
import json

# Water
gamma = 7.1  # [1]
pi_inf = 306.0e06  # [N/m2]
mu = 1.002e-03  # [kg/m/s]

# Surrounding flow
rho0 = 1000.0  # [kg/m3]
u0 = 0.0  # [m/s]
pres0 = 101325.0  # [N/m2]
c0 = math.sqrt(gamma * (pres0 + pi_inf) / rho0)

# von Karman-Pao spectrum from Saad et al. (2017, AIAA)
hit_nmode = 5000  # [1]
hit_alpha = 1.453  # [1]
hit_kappa_e = 40*math.sqrt(5/12)  # [1/m]
hit_uprime = 0.25  # [m/s]
hit_lscale = 0.746834/hit_kappa_e  # [m]
hit_epsilon = hit_uprime**3/hit_lscale  # [m2/s3]
hit_nu = mu / rho0  # [m2/s]
hit_kappa_eta = hit_epsilon**(1/4)*hit_nu**(-3/4)  # [1/m]
hit_domain = 0.09*2*math.pi  # [m]

# Domain size
Lx = hit_domain  # [m]
Ly = hit_domain  # [m]
Lz = hit_domain  # [m]

# Number of grid cells
Nx = 63
Ny = 63
Nz = 63

# Grid spacing
dx = Lx / float(Nx)  # [m]
dy = Ly / float(Ny)  # [m]
dz = Lz / float(Nz)  # [m]

#
hit_kappa_max = math.pi/dx
hit_kappa_min = 2.0*math.pi/hit_domain

# Time advancement
cfl = 0.5
T = 1.0  # [s]
dt = cfl * dx / (u0 + c0)  # [s]
Ntfinal = int(T / dt)
Ntstart = 0
Nfiles = 30
t_save = int(math.ceil((Ntfinal - 0) / float(Nfiles)))
Nt = t_save * Nfiles
t_step_start = Ntstart
t_step_stop = int(Nt)

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": Lx,
            "y_domain%beg": 0.0,
            "y_domain%end": Ly,
            "z_domain%beg": 0.0,
            "z_domain%end": Lz,
            "m": Nx,
            "n": Ny,
            "p": Nz,
            "dt": dt,
            "t_step_start": t_step_start,
            "t_step_stop": t_step_stop,
            "t_step_save": t_save,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "num_fluids": 1,
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-40,
            "weno_Re_flux": "F",
            "wenoz": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            "bc_z%beg": -1,
            "bc_z%end": -1,
            "viscous": "T",
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "cons_vars_wrt": "F",
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            "fd_order": 1,
            "omega_wrt(1)": "T",
            "omega_wrt(2)": "T",
            "omega_wrt(3)": "T",
            "qm_wrt": "T",
            "liutex_wrt": "T",
            # Patch 1
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": Lx / 2.0,
            "patch_icpp(1)%y_centroid": Ly / 2.0,
            "patch_icpp(1)%z_centroid": Lz / 2.0,
            "patch_icpp(1)%length_x": Lx,
            "patch_icpp(1)%length_y": Ly,
            "patch_icpp(1)%length_z": Lz,
            "patch_icpp(1)%alpha_rho(1)": rho0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": pres0,
            # HIT perturbation
            "hit_perturb": "T",
            "hit_nmode": hit_nmode,
            "hit_alpha": hit_alpha,
            "hit_kappa_e": hit_kappa_e,
            "hit_uprime": hit_uprime,
            "hit_kappa_eta": hit_kappa_eta,
            "hit_kappa_min": hit_kappa_min,
            "hit_kappa_max": hit_kappa_max,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%pi_inf": gamma * pi_inf / (gamma - 1.0),
            "fluid_pp(1)%Re(1)": 1.0 / mu,
        }
    )
)
