"""Benchmark driver (prompts/001 step 8)."""
from .constants import C_BASELINE, IRREPS


def run(irreps=IRREPS, energies=(0.2700, 0.3100, 0.3599), repeats=5, warmup=1, threads=1,
        backends=("oracle", "fast"), c=C_BASELINE, out_json=None):
    raise NotImplementedError("qc3.bench.run: prompts/001 step 8")


def main(argv=None):
    return run()


if __name__ == "__main__":
    main()
