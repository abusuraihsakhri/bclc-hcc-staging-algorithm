#!/usr/bin/env python3
"""Longitudinal BCLC tracking and OPTN HCC downstaging screening.

The OPTN helper is a screening aid only. Current transplant policy and transplant-
center review remain authoritative.
"""
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

from bclc_staging import stage_bclc

STAGE_ORDER = ["0", "A", "B", "C", "D"]
AFP_HIGH_THRESHOLD = 1000.0
AFP_POST_TREATMENT_THRESHOLD = 500.0


@dataclass
class HCCAssessment:
    date: str
    tumor_size_cm: float
    tumor_count: int = 1
    total_tumor_diameter_cm: Optional[float] = None
    vascular_invasion: bool = False
    extrahepatic_spread: bool = False
    ecog_ps: int = 0
    ecog_cancer_related: bool = True
    child_pugh: str = "A"
    liver_decompensation: bool = False
    transplant_candidate: Optional[bool] = None
    afp_ng_ml: float = 0.0
    treatment_started: str = ""

    def __post_init__(self):
        if self.total_tumor_diameter_cm is None:
            self.total_tumor_diameter_cm = self.tumor_size_cm * max(self.tumor_count, 1)


def bclc_stage(assessment: HCCAssessment) -> str:
    """Use the repository's single authoritative BCLC implementation."""
    return stage_bclc(
        tumor_count=assessment.tumor_count,
        tumor_size_cm=assessment.tumor_size_cm,
        child_pugh_class=assessment.child_pugh,
        ecog_ps=assessment.ecog_ps,
        ecog_cancer_related=assessment.ecog_cancer_related,
        vascular_invasion=assessment.vascular_invasion,
        extrahepatic_spread=assessment.extrahepatic_spread,
        liver_decompensation=assessment.liver_decompensation,
        transplant_candidate=assessment.transplant_candidate,
    )["bclc_stage"]


def milan(assessment: HCCAssessment) -> bool:
    if assessment.vascular_invasion or assessment.extrahepatic_spread:
        return False
    if assessment.tumor_count == 1:
        return 0 < assessment.tumor_size_cm <= 5.0
    return 2 <= assessment.tumor_count <= 3 and 0 < assessment.tumor_size_cm <= 3.0


def optn_downstaging_screen(assessment: HCCAssessment) -> Dict[str, Any]:
    """Screen lesion burden against OPTN Policy 9.5.I.iii (effective 2025-12-10)."""
    count = assessment.tumor_count
    largest = assessment.tumor_size_cm
    total = float(assessment.total_tumor_diameter_cm or 0.0)

    if count == 1:
        size_ok = 5.0 < largest <= 8.0
        detail = "one lesion >5 cm and <=8 cm"
    elif count in {2, 3}:
        size_ok = 3.0 < largest <= 5.0 and total <= 8.0
        detail = "2-3 lesions, at least one >3 cm, each <=5 cm, total diameter <=8 cm"
    elif count in {4, 5}:
        size_ok = largest < 3.0 and total <= 8.0
        detail = "4-5 lesions, each <3 cm, total diameter <=8 cm"
    else:
        size_ok = False
        detail = "lesion count outside standardized downstaging criteria"

    no_spread = not assessment.vascular_invasion and not assessment.extrahepatic_spread
    afp_le_1000 = assessment.afp_ng_ml <= AFP_HIGH_THRESHOLD

    return {
        "optn_downstaging_tumor_burden_match": size_ok,
        "no_macrovascular_invasion_or_extrahepatic_spread": no_spread,
        "afp_le_1000": afp_le_1000,
        "screen_positive": size_ok and no_spread and afp_le_1000,
        "detail": detail,
        "policy_note": (
            "OPTN policy also requires post-treatment imaging and T2-stage criteria. "
            "AFP >1000 ng/mL has a separate pathway requiring reduction below 500 ng/mL."
        ),
    }


# Backward-compatible function name.
unos_downstaging_eligible = optn_downstaging_screen


@dataclass
class StagingEpisode:
    assessment: HCCAssessment
    stage: str
    milan_ok: bool


def track_stages(assessments: List[HCCAssessment]) -> Dict[str, Any]:
    episodes = [
        StagingEpisode(item, bclc_stage(item), milan(item))
        for item in assessments
    ]
    transitions = []
    alerts = []

    for previous, current in zip(episodes, episodes[1:]):
        prev_idx = STAGE_ORDER.index(previous.stage)
        curr_idx = STAGE_ORDER.index(current.stage)
        if curr_idx > prev_idx:
            kind = "progression"
            alerts.append({
                "severity": "HIGH" if current.stage in {"C", "D"} else "WARNING",
                "title": f"BCLC progression {previous.stage}->{current.stage}",
                "date": current.assessment.date,
                "recommendation": "Multidisciplinary reassessment of imaging, liver status, and treatment options",
            })
        elif curr_idx < prev_idx:
            kind = "lower_stage_on_reassessment"
            alerts.append({
                "severity": "INFO",
                "title": f"Lower BCLC stage {previous.stage}->{current.stage}",
                "date": current.assessment.date,
                "recommendation": "Confirm response versus measurement or input differences",
            })
        else:
            kind = "stable"

        transitions.append({
            "from_date": previous.assessment.date,
            "to_date": current.assessment.date,
            "from_stage": previous.stage,
            "to_stage": current.stage,
            "kind": kind,
            "treatment_between": current.assessment.treatment_started or "none documented",
        })

    return {
        "tool": "bclc-hcc-staging-algorithm/longitudinal",
        "bclc_version": "2026",
        "timeline": [
            {
                "date": episode.assessment.date,
                "stage": episode.stage,
                "milan": episode.milan_ok,
                "afp_ng_ml": episode.assessment.afp_ng_ml,
                "ecog_ps": episode.assessment.ecog_ps,
            }
            for episode in episodes
        ],
        "transitions": transitions,
        "alerts": alerts,
        "current_downstaging_screen": (
            optn_downstaging_screen(assessments[-1])
            if assessments and not milan(assessments[-1])
            else None
        ),
    }
