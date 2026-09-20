#!/usr/bin/env python3
"""Dependency-free stress/smoke simulator for the BCLC staging engine."""
import random
import sys
import time

from bclc_staging import stage_bclc


def run_simulation(iterations: int = 100, seed: int = 42) -> dict:
    if iterations < 1:
        raise ValueError("iterations must be >= 1")

    rng = random.Random(seed)
    counts = {stage: 0 for stage in ["0", "A", "B", "C", "D"]}
    validation_blocks = 0
    start = time.perf_counter()

    for index in range(iterations):
        decomp = rng.random() < 0.12
        case = {
            "tumor_count": rng.randint(1, 7),
            "tumor_size_cm": round(rng.uniform(0.5, 10.0), 2),
            "child_pugh_class": rng.choice(["A", "A", "B", "C"]),
            "ecog_ps": rng.randint(0, 4),
            "portal_vein_invasion": rng.random() < 0.15,
            "extrahepatic_spread": rng.random() < 0.12,
            "liver_decompensation": decomp,
            "transplant_candidate": rng.choice([True, False]) if decomp else None,
        }
        result = stage_bclc(**case)
        counts[result["bclc_stage"]] += 1

        if (index + 1) % 25 == 0:
            try:
                stage_bclc(tumor_count=0, tumor_size_cm=1.0)
            except ValueError:
                validation_blocks += 1

    elapsed = time.perf_counter() - start
    summary = {
        "iterations": iterations,
        "seed": seed,
        "elapsed_seconds": round(elapsed, 6),
        "cases_per_second": round(iterations / max(elapsed, 1e-9), 1),
        "stage_counts": counts,
        "invalid_input_checks_blocked": validation_blocks,
    }
    print(summary)
    return summary


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    run_simulation(n)
