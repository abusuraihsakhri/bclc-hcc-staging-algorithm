#!/usr/bin/env python3
"""
BCLC HCC Longitudinal Stage Tracker & Transplant Downstaging Engine
Serial BCLC staging with transition classification (progression/stable/downstage),
UNOS-DS downstaging eligibility, and alert escalation rules.

Zero-dependency. Author: Dr. Abu Suraih Sakhri. License: MIT.
"""
import argparse
import json
import sys
from dataclasses import dataclass, field
from datetime import datetime
from typing import Dict, Any, List, Optional

STAGE_ORDER = ["0", "A", "B", "C", "D"]

# UNOS downstaging protocol (UCSF/UNOS-DS): AFP <= 1000 required
UNOS_DS_RULES = {
    "1_tumor": {"min_cm": 5.0, "max_cm": 8.0},
    "2_3_tumors": {"largest_max_cm": 5.0, "any_min_cm": 3.0, "total_max_cm": 8.0},
    "4_5_tumors": {"each_max_cm": 3.0, "total_max_cm": 8.0},
}
AFP_DOWNSTAGING_MAX = 1000.0   # ng/mL at presentation
AFP_POST_TARGET = 400.0        # must fall <=400 to list


@dataclass
class HCCAssessment:
    """One timepoint of HCC tumor burden + liver function."""
    date: str
    tumor_size_cm: float = 0.0       # largest nodule
    tumor_count: int = 1
    total_tumor_diameter_cm: float = None  # sum of diameters; defaults to size*count
    vascular_invasion: bool = False
    extrahepatic_spread: bool = False
    ecog_ps: int = 0
    child_pugh: str = "A"
    afp_ng_ml: float = 0.0
    treatment_started: str = ""

    def __post_init__(self):
        if self.total_tumor_diameter_cm is None:
            self.total_tumor_diameter_cm = self.tumor_size_cm * max(1, self.tumor_count)


def bclc_stage(a: HCCAssessment) -> str:
    """EASL BCLC decision tree."""
    if a.child_pugh == "C" or a.ecog_ps >= 2:
        return "D"
    if a.extrahepatic_spread or a.vascular_invasion or a.ecog_ps == 1:
        return "C"
    if a.child_pugh == "B":
        if a.tumor_count <= 3 and a.tumor_size_cm <= 3.0 and a.ecog_ps == 0:
            return "B"
        return "B"
    # Child A, ECOG 0
    if a.tumor_count == 1 and a.tumor_size_cm <= 2.0:
        return "0"
    if a.tumor_count <= 3 and a.tumor_size_cm <= 3.0:
        return "A"
    if a.tumor_count == 1 and a.tumor_size_cm <= 5.0:
        return "A"
    return "B"


def milan(a: HCCAssessment) -> bool:
    if a.tumor_count == 1:
        return a.tumor_size_cm <= 5.0
    return 2 <= a.tumor_count <= 3 and a.tumor_size_cm <= 3.0


def unos_downstaging_eligible(a: HCCAssessment) -> Dict[str, Any]:
    """UNOS-DS criteria for tumors beyond Milan but within downstaging limits."""
    reasons = []
    ok_size = False
    if a.tumor_count == 1:
        r = UNOS_DS_RULES["1_tumor"]
        ok_size = r["min_cm"] <= a.tumor_size_cm <= r["max_cm"]
        reasons.append(f"single {a.tumor_size_cm:.1f}cm in [5-8]: {'PASS' if ok_size else 'FAIL'}")
    elif a.tumor_count in (2, 3):
        r = UNOS_DS_RULES["2_3_tumors"]
        ok_size = (a.tumor_size_cm >= r["any_min_cm"] and a.tumor_size_cm <= r["largest_max_cm"]
                   and a.total_tumor_diameter_cm <= r["total_max_cm"])
        reasons.append(f"2-3 nodules largest<=5 total<=8 with one>=3: {'PASS' if ok_size else 'FAIL'}")
    elif a.tumor_count in (4, 5):
        r = UNOS_DS_RULES["4_5_tumors"]
        ok_size = a.tumor_size_cm <= r["each_max_cm"] and a.total_tumor_diameter_cm <= r["total_max_cm"]
        reasons.append(f"4-5 nodules each<=3 total<=8: {'PASS' if ok_size else 'FAIL'}")
    else:
        reasons.append("more than 5 nodules: outside UNOS-DS")

    afp_ok = a.afp_ng_ml <= AFP_DOWNSTAGING_MAX
    no_vascular = not a.vascular_invasion and not a.extrahepatic_spread
    eligible = ok_size and afp_ok and no_vascular
    return {
        "unos_ds_eligible": eligible,
        "checks": {
            "size_within_protocol": ok_size,
            "afp_le_1000": afp_ok,
            "no_vascular_or_mets": no_vascular,
        },
        "detail": reasons,
        "next_step": ("Locoregional downstaging (TACE/TARE), re-image q8 weeks; "
                      "list when within Milan AND AFP<=400")
                      if eligible else "Not a downstaging candidate",
    }


@dataclass
class StagingEpisode:
    assessment: HCCAssessment
    stage: str
    milan_ok: bool


def track_stages(assessments: List[HCCAssessment]) -> Dict[str, Any]:
    episodes = [StagingEpisode(a, bclc_stage(a), milan(a)) for a in assessments]
    transitions = []
    alerts = []
    for prev, curr in zip(episodes, episodes[1:]):
        d_prev, d_curr = STAGE_ORDER.index(prev.stage), STAGE_ORDER.index(curr.stage)
        if d_curr > d_prev:
            kind = "progression"
            alerts.append({
                "severity": "CRITICAL_ACTION_REQUIRED" if curr.stage in ("C", "D") else "WARNING",
                "title": f"BCLC progression {prev.stage}->{curr.stage} on {curr.assessment.date}",
                "recommendation": _transition_action(prev.stage, curr.stage),
            })
        elif d_curr < d_prev:
            kind = "downstaged"
            alerts.append({
                "severity": "INFO",
                "title": f"Response-induced downstaging {prev.stage}->{curr.stage}",
                "recommendation": "Re-assess transplant candidacy / curative options",
            })
        else:
            kind = "stable"

        ps_drop = prev.assessment.ecog_ps < curr.assessment.ecog_ps
        if ps_drop and d_curr - d_prev >= 1:
            alerts.append({
                "severity": "WARNING",
                "title": f"ECOG decline {prev.assessment.ecog_ps}->{curr.assessment.ecog_ps} changed stage",
                "recommendation": "Re-baseline performance status; verify imaging before therapy switch",
            })
        if prev.milan_ok and not curr.milan_ok and curr.stage != "D":
            alerts.append({
                "severity": "CRITICAL_ACTION_REQUIRED",
                "title": f"Milan eligibility threatened on {curr.assessment.date}",
                "recommendation": "Urgent transplant committee review; consider downstaging bridge",
            })
        if curr.stage == "D":
            alerts.append({
                "severity": "ADVISORY",
                "title": "BCLC-D reached",
                "recommendation": "Palliative care consultation trigger",
            })

        transitions.append({
            "from_date": prev.assessment.date, "to_date": curr.assessment.date,
            "from_stage": prev.stage, "to_stage": curr.stage,
            "kind": kind, "treatment_between": curr.assessment.treatment_started or "none documented",
        })

    return {
        "tool": "bclc-hcc-staging-algorithm/enriched",
        "timeline": [
            {"date": e.assessment.date, "stage": e.stage, "milan": e.milan_ok,
             "afp": e.assessment.afp_ng_ml, "ecog": e.assessment.ecog_ps}
            for e in episodes
        ],
        "transitions": transitions,
        "alerts": alerts,
        "current_downstaging_eval": (
            unos_downstaging_eligible(assessments[-1])
            if assessments and not milan(assessments[-1]) else None
        ),
    }


def _transition_action(frm: str, to: str) -> str:
    table = {
        ("A", "B"): "Progression beyond Milan window: escalate TACE planning, re-check AFP velocity",
        ("0", "A"): "Growth within early stage: confirm with LI-RADS 5 imaging, plan ablation/resection",
        ("B", "C"): "Portal invasion or PS decline: switch to systemic therapy (atezo/bev)",
        ("C", "D"): "Liver decompensation or PS>=2: withdraw systemic therapy, hospice evaluation",
    }
    return table.get((frm, to), "Tumor board review of stage transition")


if __name__ == "__main__":
    timeline = [
        HCCAssessment(date="2026-01-10", tumor_size_cm=2.8, tumor_count=2, ecog_ps=0,
                      child_pugh="A", afp_ng_ml=45, treatment_started=""),
        HCCAssessment(date="2026-04-15", tumor_size_cm=4.9, tumor_count=3, ecog_ps=0,
                      child_pugh="A", afp_ng_ml=210, treatment_started="TACE cycle 1"),
        HCCAssessment(date="2026-08-01", tumor_size_cm=6.2, tumor_count=3, ecog_ps=1,
                      child_pugh="B", afp_ng_ml=890, vascular_invasion=False,
                      extrahepatic_spread=False, treatment_started="TACE cycle 2"),
    ]
    print(json.dumps(track_stages(timeline), indent=2))
