"""Decision D7 (Digonto, 2026-10-22; prompts/002 sub-step 7.8): the production p-wave K2^{-1} form is "q0_consistent"
(HadSpec-consistent, docs/METHODS.md section 5), recorded in every result's metadata.  "coded_v2" stays selectable.

This test builds its Sector directly from qc3 (NOT through m1_helpers.sector, which pins "coded_v2" for the locked
comparisons against coded_v2 anchors and labels).
"""
import numpy as np

from m1_helpers import C_ZERO, LBYAS, anchor_points


def test_production_default_k2_pwave_form_is_q0_consistent():
    from qc3 import constants, fast, oracle
    from qc3.sector import Sector

    s = Sector.production("100_A2", Lbyas=LBYAS)
    assert s.k2_pwave_form == "q0_consistent"
    assert Sector("110_A2").k2_pwave_form == "q0_consistent"          # the dataclass default itself
    assert getattr(constants, "K2_PWAVE_FORM_DEFAULT") == "q0_consistent"
    assert set(oracle.K2_PWAVE_FORMS) == {"coded_v2", "q0_consistent"}  # both remain selectable
    assert Sector.production("100_A2", k2_pwave_form="coded_v2").k2_pwave_form == "coded_v2"

    p = anchor_points(small_only=True)[0]
    r = fast.evaluate_point(s if p["irrep"] == "100_A2" else Sector.production(p["irrep"], Lbyas=LBYAS), p["Ecm"], C_ZERO,
                            want=("K2inv", "Q"))
    assert r.metadata["k2_pwave_form"] == "q0_consistent"
    # the default really is the other form: the channel-1 p-wave diagonal of K2^{-1} differs from coded_v2
    r2 = fast.evaluate_point(Sector.production(p["irrep"], Lbyas=LBYAS, k2_pwave_form="coded_v2"), p["Ecm"], C_ZERO,
                             want=("K2inv", "Q"))
    assert r2.metadata["k2_pwave_form"] == "coded_v2"
    d1 = np.real(np.diag(np.asarray(r.K2inv)))[:r.N1]
    d2 = np.real(np.diag(np.asarray(r2.K2inv)))[:r2.N1]
    pwave = np.asarray(r.basis["ell"])[:r.N1] == 1
    assert np.array_equal(d1[~pwave], d2[~pwave])        # s-wave entries identical
    assert np.any(d1[pwave] != d2[pwave])                # p-wave entries differ
