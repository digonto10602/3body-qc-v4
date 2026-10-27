import sys; sys.path.insert(0,'tests/acceptance/m1')
import numpy as np
from m1_helpers import anchor_points, C_ZERO, oracle_point
p=[q for q in anchor_points() if q['irrep']=='200_A2' and abs(q['Ecm']-0.34738)<1e-9][0]
r=oracle_point(p['irrep'],p['Ecm'],C_ZERO)
print(r.lnabsdet, p['K3_zero']['lnabsdet'], r.lnabsdet-p['K3_zero']['lnabsdet'], r.metadata.get('lnabsdet_eigen_lu'))
