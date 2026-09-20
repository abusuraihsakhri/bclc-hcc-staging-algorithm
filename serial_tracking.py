"""Longitudinal BCLC stage tracking."""
from dataclasses import dataclass
from typing import List


STAGE_ORDER = {"0": 0, "A": 1, "B": 2, "C": 3, "D": 4}


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
        stage = str(bclc_stage).upper()
        if stage not in STAGE_ORDER:
            raise ValueError(f"Invalid BCLC stage: {bclc_stage!r}")
        rec = BCLCStagingRecord(
            patient_id=patient_id,
            date=date,
            bclc_stage=stage,
            **kwargs,
        )
        if not 0 <= rec.ecog_ps <= 4:
            raise ValueError("ECOG PS must be 0-4")
        self.records.append(rec)
        return rec

    def get_timeline(self, patient_id: str) -> List[BCLCStagingRecord]:
        return sorted(
            [record for record in self.records if record.patient_id == patient_id],
            key=lambda record: record.date,
        )

    def detect_progression(self, patient_id: str) -> List[dict]:
        timeline = self.get_timeline(patient_id)
        alerts = []
        for previous, current in zip(timeline, timeline[1:]):
            prev = STAGE_ORDER[previous.bclc_stage]
            curr = STAGE_ORDER[current.bclc_stage]
            if curr == prev:
                continue
            alerts.append({
                "type": "BCLC_PROGRESSION" if curr > prev else "BCLC_LOWER_STAGE",
                "patient_id": patient_id,
                "from_stage": previous.bclc_stage,
                "to_stage": current.bclc_stage,
                "date": current.date,
            })
        return alerts

    def get_cohort_summary(self) -> dict:
        distribution = {}
        for record in self.records:
            distribution[record.bclc_stage] = distribution.get(record.bclc_stage, 0) + 1
        return {"total": len(self.records), "stage_distribution": distribution}

    def detect_ecog_change_alert(self, patient_id: str) -> List[dict]:
        timeline = self.get_timeline(patient_id)
        alerts = []
        for previous, current in zip(timeline, timeline[1:]):
            if current.ecog_ps <= previous.ecog_ps:
                continue
            alerts.append({
                "type": "ECOG_SEVERE_IMPAIRMENT" if current.ecog_ps >= 3 else "ECOG_DECLINE",
                "patient_id": patient_id,
                "from_ecog": previous.ecog_ps,
                "to_ecog": current.ecog_ps,
                "date": current.date,
                "note": "Confirm whether ECOG change is HCC-attributable before using it for BCLC staging.",
            })
        return alerts
