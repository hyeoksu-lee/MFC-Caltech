"""Input contracts for the fixed-checkpoint CPU/OpenACC POEM integration."""

import importlib.util
import math
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from .case import Case
from .case_validator import CaseValidator
from .test_case_validator import BASE, ConstraintTestCase

POEM = {
    **BASE,
    "poem_euler": "T",
    "poem_distribution": "T",
    "poem_fixed_density": "T",
    "bubbles_euler": "T",
    "polytropic": "T",
    "adv_n": "T",
    "nb": 1,
    "polydisperse": "F",
    "bubble_model": 2,
    "adap_dt": "T",
    "time_stepper": 3,
    "patch_icpp(1)%alpha(1)": 1e-5,
    "patch_icpp(1)%r0": 1.0,
    "patch_icpp(1)%v0": 0.0,
    "bub_pp%R0ref": 1.0,
    "bub_pp%p0ref": 1.0,
    "bub_pp%rho0ref": 1.0,
    "bub_pp%gam_g": 1.4,
    "bub_pp%ss": 0.0728 / (101325.0 * 1e-5),
    "bub_pp%pv": 2338.8 / 101325.0,
    "bub_pp%mu_l": 0.001 / (1000.0 * 1e-5 * math.sqrt(101325.0 / 1000.0)),
    "fluid_pp(1)%gamma": 1.0 / 6.1,
    "fluid_pp(1)%pi_inf": 7.1 * (3.06e8 / 101325.0) / 6.1,
    "poly_sigma": 0.7,
}


class PoemConstraints(ConstraintTestCase):
    def test_openacc_allowed_but_openmp_and_reduced_precision_rejected(self):
        for gpu, single, mixed, accepted in (("acc", False, False, True), ("mp", False, False, False), ("acc", True, False, False), ("acc", False, True, False)):
            with self.subTest(gpu=gpu, single=single, mixed=mixed):
                with patch("mfc.case_validator.CFG", return_value=SimpleNamespace(gpu=gpu, single=single, mixed=mixed)):
                    validator = CaseValidator(POEM)
                    validator.check_poem()
                    self.assertEqual(not validator.errors, accepted)

    def test_supported_case(self):
        self.assertAccepts(POEM)

    def test_no_fictitious_bins(self):
        self.assertRejects({**POEM, "nb": 32}, "not radius bins")

    def test_number_density_required(self):
        self.assertRejects({**POEM, "adv_n": "F"}, "poem_euler requires adv_n")

    def test_source_splitting_required(self):
        self.assertRejects({**POEM, "adap_dt": "F"}, "split source path")

    def test_fixed_source_physics_required(self):
        self.assertRejects({**POEM, "poem_fixed_density": "F"}, "poem_fixed_density = T")

    def test_nonreference_initialization_rejected(self):
        self.assertRejects({**POEM, "patch_icpp(1)%r0": 0.8}, "reference bubble radii")

    def test_multidimensional_case_rejected(self):
        self.assertRejects({**POEM, "n": 50}, "requires a 1D domain")

    def test_compact_diagnostics_supported(self):
        self.assertAccepts({**POEM, "poem_diagnostics": "T"})

    def test_radius_probe_rejected(self):
        self.assertRejects({**POEM, "probe_wrt": "T"}, "radius-bin probes")

    def test_radius_postprocess_rejected(self):
        validator = CaseValidator(POEM)
        validator.validate_post_process()
        self.assertIn("POEM post_process radius-bin output is unsupported", "\n".join(validator.errors))

    def test_mismatched_physics_rejected(self):
        self.assertRejects({**POEM, "bub_pp%gam_g": 1.0}, "frozen checkpoint requires explicit bub_pp%gam_g")

    def test_missing_physics_rejected(self):
        self.assertRejects({k: v for k, v in POEM.items() if k != "bub_pp%pv"}, "frozen checkpoint requires explicit bub_pp%pv")

    def test_fixed_density_direct_requires_split_source(self):
        self.assertRejects({**POEM, "poem_euler": "F", "adap_dt": "F"}, "poem_fixed_density requires the split source")

    def test_defaults_leave_direct_unchanged(self):
        self.assertAccepts({key: value for key, value in POEM.items() if not key.startswith("poem_")})


class GeneratedPoemCases(unittest.TestCase):
    """Exercise the actual coupled-case generator through the native validators."""

    def test_all_kinds_and_models_validate_for_both_targets(self):
        path = Path(__file__).resolve().parents[3] / "POEM/toolchain/workflows/validation/mfc_integration/coupled.py"
        spec = importlib.util.spec_from_file_location("poem_coupled_case_regression", path)
        cases = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cases)
        for kind in cases.KINDS:
            for model in cases.MODELS:
                for target in ("pre_process", "simulation"):
                    with self.subTest(kind=kind, model=model, target=target):
                        case = Case(cases.coupled_case(kind, model))
                        case.validate_params()
                        CaseValidator(case.params).validate(target)
