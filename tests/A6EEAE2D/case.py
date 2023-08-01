#!/usr/bin/env python3
#
# tests/A6EEAE2D/case.py:
# 2D -> Bubbles -> Monopole -> QBMM -> bubble_model=3 -> Non-polytropic

import json
import argparse

parser = argparse.ArgumentParser(
    prog="tests/A6EEAE2D/case.py",
    description="tests/A6EEAE2D/case.py: 2D -> Bubbles -> Monopole -> QBMM -> bubble_model=3 -> Non-polytropic",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
 
parser.add_argument("dict", type=str, metavar="DICT", help=argparse.SUPPRESS)

ARGS = vars(parser.parse_args())

ARGS["dict"] = json.loads(ARGS["dict"])

case = {
    "run_time_info": "T",
    "m": 49,
    "n": 39,
    "p": 0,
    "dt": 1e-06,
    "t_step_start": 0,
    "t_step_stop": 50,
    "t_step_save": 50,
    "num_patches": 3,
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "num_fluids": 1,
    "adv_alphan": "T",
    "mpp_lim": "F",
    "mixture_err": "F",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "F",
    "null_weights": "F",
    "mp_weno": "F",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "F",
    "parallel_io": "F",
    "patch_icpp(1)%pres": 1.0,
    "patch_icpp(1)%alpha_rho(1)": 0.96,
    "patch_icpp(1)%alpha(1)": 0.04,
    "patch_icpp(2)%pres": 1.0,
    "patch_icpp(2)%alpha_rho(1)": 0.96,
    "patch_icpp(2)%alpha(1)": 0.04,
    "patch_icpp(3)%pres": 1.0,
    "patch_icpp(3)%alpha_rho(1)": 0.96,
    "patch_icpp(3)%alpha(1)": 0.04,
    "fluid_pp(1)%gamma": 0.16,
    "fluid_pp(1)%pi_inf": 3515.0,
    "bubbles": "T",
    "Ca": 0.9769178386380458,
    "Web": 13.927835051546392,
    "Re_inv": 0.009954269975623245,
    "pref": 101325.0,
    "rhoref": 1000.0,
    "bubble_model": 3,
    "polytropic": "F",
    "polydisperse": "F",
    "thermal": 3,
    "R0ref": 1e-05,
    "patch_icpp(1)%r0": 1,
    "patch_icpp(1)%v0": 0,
    "patch_icpp(2)%r0": 1,
    "patch_icpp(2)%v0": 0,
    "patch_icpp(3)%r0": 1,
    "patch_icpp(3)%v0": 0,
    "qbmm": "T",
    "dist_type": 2,
    "poly_sigma": 0.3,
    "R0_type": 1,
    "sigR": 0.1,
    "sigV": 0.1,
    "rhoRV": 0.0,
    "Monopole": "T",
    "num_mono": 1,
    "Mono(1)%loc(1)": 0.5,
    "Mono(1)%mag": 1.0,
    "Mono(1)%length": 0.25,
    "Mono(1)%dir": 1.0,
    "Mono(1)%npulse": 1,
    "Mono(1)%pulse": 1,
    "cu_mpi": "F",
    "x_domain%beg": 0.0,
    "x_domain%end": 1.0,
    "y_domain%beg": 0.0,
    "y_domain%end": 1.0,
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -3,
    "bc_y%end": -3,
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%y_centroid": 0.05,
    "patch_icpp(1)%length_y": 0.1,
    "patch_icpp(2)%y_centroid": 0.45,
    "patch_icpp(2)%length_y": 0.7,
    "patch_icpp(3)%y_centroid": 0.9,
    "patch_icpp(3)%length_y": 0.2,
    "patch_icpp(1)%x_centroid": 0.5,
    "patch_icpp(1)%length_x": 1,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(2)%geometry": 3,
    "patch_icpp(2)%x_centroid": 0.5,
    "patch_icpp(2)%length_x": 1,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(3)%geometry": 3,
    "patch_icpp(3)%x_centroid": 0.5,
    "patch_icpp(3)%length_x": 1,
    "patch_icpp(3)%vel(1)": 0.0,
    "patch_icpp(3)%vel(2)": 0.0,
    "nb": 3,
    "fluid_pp(2)%gamma": 2.5,
    "fluid_pp(2)%pi_inf": 0.0,
    "fluid_pp(1)%mul0": 0.001002,
    "fluid_pp(1)%ss": 0.07275,
    "fluid_pp(1)%pv": 2338.8,
    "fluid_pp(1)%gamma_v": 1.33,
    "fluid_pp(1)%M_v": 18.02,
    "fluid_pp(1)%mu_v": 8.816e-06,
    "fluid_pp(1)%k_v": 0.019426,
    "fluid_pp(2)%gamma_v": 1.4,
    "fluid_pp(2)%M_v": 28.97,
    "fluid_pp(2)%mu_v": 1.8e-05,
    "fluid_pp(2)%k_v": 0.02556,
    "Mono(1)%loc(2)": 0.5
}
mods = {}

if "post_process" in ARGS["dict"]["targets"]:
    mods = {
        'parallel_io'  : 'T', 'cons_vars_wrt'   : 'T',
        'prim_vars_wrt': 'T', 'alpha_rho_wrt(1)': 'T',
        'rho_wrt'      : 'T', 'mom_wrt(1)'      : 'T',
        'vel_wrt(1)'   : 'T', 'E_wrt'           : 'T',
        'pres_wrt'     : 'T', 'alpha_wrt(1)'    : 'T',
        'gamma_wrt'    : 'T', 'heat_ratio_wrt'  : 'T',
        'pi_inf_wrt'   : 'T', 'pres_inf_wrt'    : 'T',
        'c_wrt'        : 'T',
    }
        
    if case['p'] != 0:
        mods['fd_order']  = 1
        mods['omega_wrt'] = 'T'

print(json.dumps({**case, **mods}))