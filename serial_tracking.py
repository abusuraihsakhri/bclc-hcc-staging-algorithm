"""Longitudinal BCLC stage tracking."""
from dataclasses import dataclass
from typing import List


@dataclass
class BCLCStagingRecord:
    patient_id: str
    date: str
    bclc_stage: str
    tumor_size_cm: float = 0.0
    tumor_count: int = 1
    ecog_ps: int = 0
    child_pugh: str = "A"
    milan_eligible: bool = False
    treatment: str = ""
    notes: str = ""


class BCLCTracker:
    def __init__(self):
        self.records: List[BCLCStagingRecord] = []

    def add_staging(self, patient_id: str, date: str, bclc_stage: str, **kwargs) -> BCLCStagingRecord:
        rec = BCLCStagingRecord(
            patient_id=patient_id, date=date, bclc_stage=bclc_stage.upper(), **kwargs,
        )
        self.records.append(rec)
        return rec

    def get_timeline(self, patient_id: str) -> List[BCLCStagingRecord]:
        return sorted(
            [r for r in self.records if r.patient_id == patient_id],
            key=lambda r: r.date,
        )

    def detect_progression(self, patient_id: str) -> List[dict]:
        timeline = self.get_timeline(patient_id)
        if len(timeline) < 2:
            return []
        stage_order = {"0": 0, "A": 1, "B": 2, "C": 3, "D": 4}
        alerts = []
        for i in range(1, len(timeline)):
            prev = stage_order.get(timeline[i - 1].bclc_stage, 0)
            curr = stage_order.get(timeline[i].bclc_stage, 0)
            if curr > prev:
                alerts.append({
                    "type": "BCLC_PROGRESSION",
                    "patient_id": patient_id,
                    "from_stage": timeline[i - 1].bclc_stage,
                    "to_stage": timeline[i].bclc_stage,
                    "date": timeline[i].date,
                })
            elif curr < prev:
                alerts.append({
                    "type": "BCLC_DOWNSTAGING",
                    "patient_id": patient_id,
                    "from_stage": timeline[i - 1].bclc_stage,
                    "to_stage": timeline[i].bclc_stage,
                    "date": timeline[i].date,
                })
        return alerts

    def get_cohort_summary(self) -> dict:
        if not self.records:
            return {"total": 0}
        stages = {}
        for r in self.records:
            stages[r.bclc_stage] = stages.get(r.bclc_stage, 0) + 1
        return {"total": len(self.records), "stage_distribution": stages}

    def detect_ecog_change_alert(self, patient_id: str) -> List[dict]:
        timeline = self.get_timeline(patient_id)
        if len(timeline) < 2:
            return []
        alerts = []
        for i in range(1, len(timeline)):
            prev_ecog = timeline[i - 1].ecog_ps
            curr_ecog = timeline[i].ecog_ps
            if curr_ecog >= 2 and prev_ecog < 2:
                alerts.append({
                    "type": "ECOG_DECLINE_TO_TERMINAL",
                    "patient_id": patient_id,
                    "from_ecog": prev_ecog,
                    "to_ecog": curr_ecog,
                    "date": timeline[i].date,
                    "action": "Consider BCLC-D reassessment and palliative care referral",
                })
            elif curr_ecog > prev_ecog:
                alerts.append({
                    "type": "ECOG_DECLINE",
                    "patient_id": patient_id,
                    "from_ecog": prev_ecog,
                    "to_ecog": curr_ecog,
                    "date": timeline[i].date,
                })
        return alerts
