# L20 full Ecm-cutoff fit suite

L24 was not resumed by this task. All runnable cases use L20 accepted sectors and projected-basis hot-window caches only.

| cutoff | levels | model_found | chi2 | dof | chi2/dof | K3iso0 | K3iso1 | K3B | K3E | status |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 0.315 | 6 | 6/6 | 3.262140074829519e-05 | 2 | 1.6310700374147595e-05 | 151422.46172441466 | -972421.14060757787 | 347174.05548116949 | -1226756.7068845264 | PASS |
| 0.325 | 9 | 9/9 | 9.530896622369882e-05 | 5 | 1.9061793244739764e-05 | 125968.22184657381 | -972421.14060757787 | 347174.05548116949 | -1226756.7068845264 | PASS |
| 0.335 | 12 | 12/12 | 0.00012832616631556945 | 8 | 1.6040770789446181e-05 | 32040.841899235049 | -1074716.6512728098 | 347174.05548116943 | -1226756.7068845264 | PASS |
| 0.345 | 14 | 14/14 | 0.00016711008405490148 | 10 | 1.6711008405490147e-05 | 40607.515418119001 | -1067567.8661876961 | 376187.1404411546 | -1226756.7068845264 | PASS |

## Interpretation

The validated accepted-window package covers through 0.335. Cutoffs 0.345 and 0.355 are reported as blocked when extra lattice levels lack accepted windows; no fake roots are inserted. Among runnable cutoffs, inspect parameter drift and covariance together before production interpretation.
