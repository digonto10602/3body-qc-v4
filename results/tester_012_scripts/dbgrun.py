import sys
sys.path.insert(0, "/tmp/claude-0/-home-user-3body-qc-v4/756827e6-de12-5b97-9946-026a972abb7d/scratchpad/dbg/pkg"); sys.path.insert(0, "/home/user/3body-qc-v4/tests/acceptance/m1")
from qc3 import fast
import qc3; assert "dbg" in qc3.__file__, qc3.__file__
from m1_helpers import sector, C_BASELINE
irrep, E, mode = sys.argv[1], float(sys.argv[2]), sys.argv[3]
fast.evaluate_point(sector(irrep), E, C_BASELINE, want=("Q",), options={"use_classes": mode == "classes"})
