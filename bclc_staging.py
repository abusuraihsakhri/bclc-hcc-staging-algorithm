#!/usr/bin/env python3
"""BCLC 2026 staging helper for hepatocellular carcinoma.

Implements stage-defining features from the BCLC 2026 strategy update for
education, research, and reproducible data processing. It is not a substitute
for multidisciplinary clinical assessment.

Reference: Reig M, et al. J Hepatol. 2026;84(3):631-654.
doi:10.1016/j.jhep.2025.10.020
"""

from __future__ import annotations

import argparse
import csv
import json
from typing import Any, Dict, Optional


BCLC_STAGES = {
    "0": {
        "name": "Very Early",
        "treatment": "Ablation or resection; transplant in selected candidates",
        "expected_survival": ">5 years after effective first treatment",
    },
    "A": {
        "name": "Early",
        "treatment": "Resection, ablation, or liver transplantation according to tumor and liver context",
        "expected_survival": ">5 years after effective first treatment",
    },
    "B": {
        "name": "Intermediate",
        "treatment": "Subgroup assessment for transplant, locoregional therapy, or systemic therapy",
        "expected_survival": ">2.5 years after effective first treatment",
    },
    "C": {
        "name": "Advanced",
        "treatment": "Systemic therapy is the usual evidence-based first option; individualize in MDT review",
        "expected_survival": "approximately 2 years after effective first treatment",
    },
    "D": {
        "name": "End Stage",
        "treatment": "Supportive/palliative care when transplantation is not an option",
        "expected_survival": "<1 year",
    },
}

MILAN_SINGLE_MAX_CM = 5.0
MILAN_MULTI_MAX_COUNT = 3
MILAN_MULTI_MAX_CM = 3.0


def stage_bclc(
    tumor_count: int = 1,
    tumor_size_cm: float = 0.0,
    child_pugh_class: str = "A",
    ecog_ps: int = 0,
    portal_vein_invasion: bool = False,
    extrahepatic_spread: bool = False,
    lymph_node_metastasis: bool = False,
    vascular_invasion: bool = False,
    liver_decompensation: bool = False,
    transplant_candidate: Optional[bool] = None,
    ecog_cancer_related: bool = True,
) -> Dict[str, Any]:
    """Assign a BCLC 2026 stage from stage-defining inputs.

    Child-Pugh class is retained for compatibility and descriptive output, but
    it does not define stage by itself. If decompensation is present,
    transplant candidacy must be supplied because BCLC-D cannot safely be
    inferred from Child-Pugh class alone.

    vascular_invasion is interpreted as macrovascular invasion.
    """
    cp = str(child_pugh_class).strip().upper()
    if cp not in {"A", "B", "C"}:
        raise ValueError("Child-Pugh class must be A, B, or C")
    if isinstance(tumor_count, bool) or not isinstance(tumor_count, int) or tumor_count < 1:
        raise ValueError("Tumor count must be an integer >= 1 for confirmed HCC")
    if (
        isinstance(tumor_size_cm, bool)
        or not isinstance(tumor_size_cm, (int, float))
        or tumor_size_cm <= 0
    ):
        raise ValueError("Largest tumor diameter must be > 0 cm for confirmed HCC")
    if isinstance(ecog_ps, bool) or not isinstance(ecog_ps, int) or not 0 <= ecog_ps <= 4:
        raise ValueError("ECOG PS must be an integer from 0 to 4")
    if transplant_candidate not in {None, True, False}:
        raise ValueError("transplant_candidate must be true, false, or null")
    if liver_decompensation and transplant_candidate is None:
        raise ValueError(
            "transplant_candidate is required when liver_decompensation=True; "
            "decompensation alone does not establish BCLC-D"
        )

    size = float(tumor_size_cm)
    macrovascular = bool(portal_vein_invasion or vascular_invasion)
    metastatic = bool(extrahepatic_spread or lymph_node_metastasis)

    stage = _determine_stage(
        count=tumor_count,
        size=size,
        ecog=ecog_ps,
        macrovascular_invasion=macrovascular,
        metastatic_spread=metastatic,
        liver_decompensation=liver_decompensation,
        transplant_candidate=transplant_candidate,
        ecog_cancer_related=ecog_cancer_related,
    )
    milan_eligible = _check_milan(
        tumor_count,
        size,
        macrovascular_invasion=macrovascular,
        extrahepatic_spread=extrahepatic_spread,
        lymph_node_metastasis=lymph_node_metastasis,
    )
    info = BCLC_STAGES[stage]

    return {
        "tool": "bclc-hcc-staging-algorithm",
        "bclc_version": "2026",
        "bclc_stage": stage,
        "stage_name": info["name"],
        "classification": f"BCLC Stage {stage} ({info['name']})",
        "treatment_allocation": info["treatment"],
        "expected_survival": info["expected_survival"],
        "median_survival_months": None,
        "five_year_survival_pct": None,
        "milan_criteria_eligible": milan_eligible,
        "transplant_eligible": transplant_candidate,
        "transplant_assessment": (
            "Externally supplied transplant candidacy"
            if transplant_candidate is not None
            else "Not determined: requires transplant-center assessment beyond Milan criteria"
        ),
        "clinical_recommendation": _recommendation(
            stage, liver_decompensation, transplant_candidate
        ),
        "inputs": {
            "tumor_count": tumor_count,
            "tumor_size_cm": size,
            "child_pugh_class": cp,
            "ecog_ps": ecog_ps,
            "ecog_cancer_related": bool(ecog_cancer_related),
            "portal_vein_invasion": bool(portal_vein_invasion),
            "vascular_invasion": bool(vascular_invasion),
            "extrahepatic_spread": bool(extrahepatic_spread),
            "lymph_node_metastasis": bool(lymph_node_metastasis),
            "liver_decompensation": bool(liver_decompensation),
            "transplant_candidate": transplant_candidate,
        },
        "limitations": [
            "Assumes confirmed HCC and current imaging-defined tumor burden.",
            "ECOG-based upstaging should reflect HCC-attributable symptoms.",
            "Treatment allocation is a stage-level summary, not a prescription.",
            "Milan criteria alone do not establish transplant eligibility.",
        ],
    }


def _determine_stage(
    *,
    count: int,
    size: float,
    ecog: int,
    macrovascular_invasion: bool,
    metastatic_spread: bool,
    liver_decompensation: bool,
    transplant_candidate: Optional[bool],
    ecog_cancer_related: bool,
) -> str:
    if ecog_cancer_related and ecog >= 3:
        return "D"
    if liver_decompensation and transplant_candidate is False:
        return "D"
    if macrovascular_invasion or metastatic_spread:
        return "C"
    if ecog_cancer_related and ecog in {1, 2}:
        return "C"
    if count == 1:
        return "0" if size <= 2.0 else "A"
    if count <= 3 and size <= 3.0:
        return "A"
    return "B"


def _check_milan(
    count: int,
    size: float,
    macrovascular_invasion: bool = False,
    extrahepatic_spread: bool = False,
    lymph_node_metastasis: bool = False,
) -> bool:
    """Check radiologic Milan tumor-burden criteria and spread exclusions."""
    if count < 1 or size <= 0:
        return False
    if macrovascular_invasion or extrahepatic_spread or lymph_node_metastasis:
        return False
    if count == 1:
        return size <= MILAN_SINGLE_MAX_CM
    return count <= MILAN_MULTI_MAX_COUNT and size <= MILAN_MULTI_MAX_CM


def _recommendation(
    stage: str,
    liver_decompensation: bool,
    transplant_candidate: Optional[bool],
) -> str:
    if liver_decompensation and transplant_candidate:
        return (
            "Decompensated liver disease with externally established transplant "
            "candidacy: prioritize transplant-center and multidisciplinary review."
        )
    if stage == "D":
        return (
            "End-stage pathway: supportive and palliative care are central when "
            "transplantation is not an option; reassess reversible contributors."
        )
    if stage == "C":
        return (
            "Advanced-stage pathway: systemic therapy is the usual evidence-based "
            "first option; regimen choice requires current guideline and liver-function review."
        )
    if stage == "B":
        return (
            "Intermediate-stage pathway: assess transplant criteria, tumor distribution, "
            "portal flow, liver reserve, and locoregional versus systemic treatment suitability."
        )
    if stage == "A":
        return (
            "Early-stage pathway: evaluate resection, ablation, and/or liver "
            "transplantation according to tumor pattern and liver status."
        )
    return (
        "Very-early-stage pathway: curative-intent ablation or resection is commonly "
        "considered; liver and transplant context may alter the preferred option."
    )


def calculate_metrics(
    tumor_size_cm: float = 0.0,
    tumor_count: int = 1,
    ecog_ps: int = 0,
    child_pugh: str = "A",
    child_pugh_class: Optional[str] = None,
    portal_invasion: bool = False,
    portal_vein_invasion: Optional[bool] = None,
    **kwargs: Any,
) -> Dict[str, Any]:
    """Backward-compatible wrapper around stage_bclc."""
    cp = child_pugh_class if child_pugh_class is not None else child_pugh
    pvi = portal_vein_invasion if portal_vein_invasion is not None else portal_invasion
    return stage_bclc(
        tumor_count=tumor_count,
        tumor_size_cm=tumor_size_cm,
        child_pugh_class=cp,
        ecog_ps=ecog_ps,
        portal_vein_invasion=pvi,
        **kwargs,
    )


def _parse_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        raise ValueError("Boolean value cannot be null")
    normalized = str(value).strip().lower()
    if normalized in {"true", "1", "yes", "y", "t"}:
        return True
    if normalized in {"false", "0", "no", "n", "f"}:
        return False
    raise ValueError(f"Invalid boolean value: {value!r}")


def _parse_optional_bool(value: Any) -> Optional[bool]:
    if value is None or str(value).strip() == "":
        return None
    return _parse_bool(value)


def _required(row: Dict[str, str], field: str, row_number: int) -> str:
    value = row.get(field)
    if value is None or str(value).strip() == "":
        raise ValueError(f"row {row_number}: missing required field '{field}'")
    return str(value).strip()


def process_batch(input_csv: str, output_csv: str) -> int:
    """Process a CSV cohort and write BCLC results plus row-level errors."""
    with open(input_csv, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)

    result_fields = [
        "bclc_stage",
        "stage_name",
        "treatment_allocation",
        "expected_survival",
        "milan_criteria_eligible",
        "transplant_eligible",
        "clinical_recommendation",
        "error",
    ]
    out_fields = list(dict.fromkeys(fieldnames + result_fields))
    out_rows = []

    for row_number, row in enumerate(rows, start=2):
        out = dict(row)
        try:
            decomp = _parse_bool(row.get("liver_decompensation", "false"))
            transplant = _parse_optional_bool(row.get("transplant_candidate", ""))
            result = stage_bclc(
                tumor_count=int(_required(row, "tumor_count", row_number)),
                tumor_size_cm=float(_required(row, "tumor_size_cm", row_number)),
                child_pugh_class=row.get("child_pugh_class", "A") or "A",
                ecog_ps=int(row.get("ecog_ps", 0) or 0),
                ecog_cancer_related=_parse_bool(row.get("ecog_cancer_related", "true")),
                portal_vein_invasion=_parse_bool(row.get("portal_vein_invasion", "false")),
                extrahepatic_spread=_parse_bool(row.get("extrahepatic_spread", "false")),
                lymph_node_metastasis=_parse_bool(row.get("lymph_node_metastasis", "false")),
                vascular_invasion=_parse_bool(row.get("vascular_invasion", "false")),
                liver_decompensation=decomp,
                transplant_candidate=transplant,
            )
            for key in result_fields[:-1]:
                out[key] = result[key]
            out["error"] = ""
        except (ValueError, TypeError, KeyError) as exc:
            for key in result_fields[:-1]:
                out[key] = ""
            out["error"] = str(exc)
        out_rows.append(out)

    with open(output_csv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=out_fields)
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"Processed {len(out_rows)} records -> {output_csv}")
    return len(out_rows)


def main(argv: Optional[list[str]] = None) -> None:
    parser = argparse.ArgumentParser(
        description="BCLC 2026 hepatocellular carcinoma staging helper"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    single = subparsers.add_parser("single", help="Stage a single confirmed HCC case")
    single.add_argument("--tumor-count", type=int, default=1)
    single.add_argument("--tumor-size-cm", type=float, required=True)
    single.add_argument("--child-pugh-class", default="A", choices=["A", "B", "C"])
    single.add_argument("--ecog-ps", type=int, default=0)
    single.add_argument("--ecog-not-cancer-related", action="store_true")
    single.add_argument("--portal-vein-invasion", action="store_true")
    single.add_argument("--extrahepatic-spread", action="store_true")
    single.add_argument("--lymph-node-metastasis", action="store_true")
    single.add_argument("--vascular-invasion", action="store_true")
    single.add_argument("--liver-decompensation", action="store_true")
    transplant = single.add_mutually_exclusive_group()
    transplant.add_argument(
        "--transplant-candidate", dest="transplant_candidate", action="store_true"
    )
    transplant.add_argument(
        "--not-transplant-candidate", dest="transplant_candidate", action="store_false"
    )
    single.set_defaults(transplant_candidate=None)

    batch = subparsers.add_parser("batch", help="Batch-process a CSV file")
    batch.add_argument("-i", "--input", required=True)
    batch.add_argument("-o", "--output", default="results.csv")

    args = parser.parse_args(argv)
    if args.command == "batch":
        process_batch(args.input, args.output)
        return

    try:
        result = stage_bclc(
            tumor_count=args.tumor_count,
            tumor_size_cm=args.tumor_size_cm,
            child_pugh_class=args.child_pugh_class,
            ecog_ps=args.ecog_ps,
            ecog_cancer_related=not args.ecog_not_cancer_related,
            portal_vein_invasion=args.portal_vein_invasion,
            extrahepatic_spread=args.extrahepatic_spread,
            lymph_node_metastasis=args.lymph_node_metastasis,
            vascular_invasion=args.vascular_invasion,
            liver_decompensation=args.liver_decompensation,
            transplant_candidate=args.transplant_candidate,
        )
    except ValueError as exc:
        parser.error(str(exc))
        return
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
