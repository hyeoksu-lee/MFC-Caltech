#!/usr/bin/env python
import math
import json

# FLUID PROPERTIES
# Water
n_tait = 7.1
B_tait = 306.0e06
rho0 = 1.0e03
mul0 = 1.002e-03
ss = 0.07275
pv = 2.3388e03

# Vapor
gamma_v = 1.33
M_v = 18.02
mu_v = 0.8816e-05
k_v = 0.019426

# Air
gamma_n = 1.4
M_n = 28.97
mu_n = 1.8e-05
k_n = 0.02556

# REFERENCE VALUES
R0ref = 50.0e-06
p0eq = 8236.0
ub0 = math.sqrt(p0eq/rho0)
x0 = 4*R0ref
p0 = 3*p0eq
u0 = math.sqrt(p0/rho0)

# BUBBLES
vf0 = 1e-5
nb = 1

# External pressure
pext = 1000*p0eq/p0

# DOMAIN
Nx = 30
Ldomain = 20.0e-03
L = Ldomain / x0
dx = L / float(Nx + 1)

# TIME STEPS
Tfinal = 0.05*(R0ref/ub0)/(x0/u0)
Nt = int(5e2 + 1)
t_save = 1
dt = Tfinal / (Nt - 1)

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": -0.5 * L,
            "x_domain%end": 0.5 * L,
            "stretch_x": "F",
            "cyl_coord": "F",
            "m": Nx,
            "n": 0,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": t_save,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            "fd_order": 1,
            # Patch 1 _ Background
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%length_x": L,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%pres": pext,
            "patch_icpp(1)%alpha_rho(1)": (1.0 - vf0),
            "patch_icpp(1)%alpha(1)": vf0,
            "patch_icpp(1)%r0": 1.0,
            "patch_icpp(1)%v0": 0.0,
            # Reference scales
            "bub_refs%rho0": rho0,
            "bub_refs%x0": x0,
            "bub_refs%u0": u0,
            "bub_refs%p0": p0,
            "bub_refs%R0ref": R0ref,
            "bub_refs%ub0": ub0,
            "bub_refs%p0eq": p0eq,
            # Bubbles
            "bubbles_euler": "T",
            "bubble_model": 2,
            "bub_ss": "T",
            "bub_visc": "T",
            # adv_n
            "adv_n": "T",
            # adap_dt
            "adap_dt": "T",
            "adap_dt_max_iters": 200,
            # Gas compression model
            "polytropic": "T",
            "thermal": 1,
            # Polydispersity
            "polydisperse": "F",
            "nb": nb,
            # QBMM
            "qbmm": "F",
            # Fluids Physical Parameters
            # Surrounding liquid
            "fluid_pp(1)%gamma": 1.0e00 / (n_tait - 1.0e00),
            "fluid_pp(1)%pi_inf": n_tait * (B_tait / p0) / (n_tait - 1.0),
            "fluid_pp(1)%ss": ss,
            "fluid_pp(1)%pv": pv,
            "fluid_pp(1)%mul0": mul0,
            # Last fluid_pp is always reserved for bubble gas state
            # if applicable
            "fluid_pp(2)%gamma": 1.0 / (gamma_n - 1.0),
            "fluid_pp(2)%pi_inf": 0.0e00,
        }
    )
)
