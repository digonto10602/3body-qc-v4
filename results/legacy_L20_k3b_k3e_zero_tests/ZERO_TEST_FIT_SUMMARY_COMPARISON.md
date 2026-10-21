# Zero-test fit comparison

| test | cutoff | levels | floating | fixed | model_found | chi2/dof | K3iso0 | K3iso1 | K3B | K3E | status |
|---|---:|---:|---|---|---|---:|---:|---:|---:|---:|---|
| fixed_K3B0_K3E0 | 0.315 | 6 | K3iso0,K3iso1 | K3B,K3E | 6/6 | 1.0542666800363224e-05 | 87875.45828284955 | -972421.1406075779 | 0.0 | 0.0 | PASS |
| fixed_K3B0_K3E0 | 0.325 | 9 | K3iso0,K3iso1 | K3B,K3E | 9/9 | 1.247421508822199e-05 | 102015.39501289788 | -972421.1406075779 | 0.0 | 0.0 | PASS |
| fixed_K3B0_K3E0 | 0.335 | 12 | K3iso0,K3iso1 | K3B,K3E | 11/12 | nan | nan | nan | 0.0 | 0.0 | BLOCKED_MODEL_FOUND |
| fixed_K3B0_K3E0 | 0.345 | 14 | K3iso0,K3iso1 | K3B,K3E | 0/14 | nan | nan | nan | 0.0 | 0.0 | BLOCKED_MODEL_FOUND |
| float_from_K3B0_K3E0 | 0.315 | 6 | K3iso0,K3iso1,K3B,K3E | none | 6/6 | 2.039453905725892e-05 | 99189.40229967263 | -972421.1406075779 | 0.0 | 0.0 | PASS |
| float_from_K3B0_K3E0 | 0.325 | 9 | K3iso0,K3iso1,K3B,K3E | none | 9/9 | 1.692892130183815e-05 | 124643.51944580376 | -972421.1406075779 | 0.0 | 0.0 | PASS |
| float_from_K3B0_K3E0 | 0.335 | 12 | K3iso0,K3iso1,K3B,K3E | none | 11/12 | nan | nan | nan | nan | nan | BLOCKED_MODEL_FOUND |
| float_from_K3B0_K3E0 | 0.345 | 14 | K3iso0,K3iso1,K3B,K3E | none | 0/14 | nan | nan | nan | nan | nan | BLOCKED_MODEL_FOUND |
