"""Production constants (CLAUDE.md invariant 3; prompts/001 section 5.1)."""

PRODUCTION: dict = {
    "xi": 3.444, "atmK": 0.09698, "atmpi": 0.06906, "eta_1": 1.0, "eta_2": 0.5,
    "alpha": 0.5, "epsilon_H": 0.0, "max_shell_num": 20, "waves_1": (0, 1), "waves_2": (0,),
    "parity": -1, "eig_tol": 0.05,
    "scatter": {"a0_Kpi": 4.04, "a1_Kpi": -43.2, "a0_KK": 4.12},
    "Q0norm": True, "z_min": 0.02,
}
# label-generation vector (research 004 section 3.6)
C_BASELINE: tuple = (73735.840894011912, -972421.14060757787, 347174.05548116949, -1226756.7068845264)
IRREPS: tuple = ("000_A1m", "100_A2", "110_A2", "111_A2", "200_A2")
NNP: dict = {"000_A1m": (0, 0, 0), "100_A2": (0, 0, 1), "110_A2": (1, 1, 0),
             "111_A2": (1, 1, 1), "200_A2": (0, 0, 2)}
GRID_E0, GRID_DE, GRID_N = 0.26310, 0.0969 / 19999, 20000
