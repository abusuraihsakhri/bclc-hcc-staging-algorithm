import pytest

from alert_escalation import BCLCAlertEngine
from serial_tracking import BCLCTracker
from simulator import run_simulation
from stage_transition_tracker import HCCAssessment, bclc_stage, optn_downstaging_screen
from stratification import stratify_treatment


def test_longitudinal_tracker_validates_stage_and_progression():
    tracker = BCLCTracker()
    tracker.add_staging("P1", "2026-01-01", "A")
    tracker.add_staging("P1", "2026-06-01", "B")
    assert tracker.detect_progression("P1")[0]["type"] == "BCLC_PROGRESSION"
    with pytest.raises(ValueError):
        tracker.add_staging("P1", "2026-07-01", "X")


def test_ecog_alert_no_longer_calls_ecog2_terminal():
    alerts = BCLCAlertEngine().evaluate(
        "P1", "C", False, 2, "A", 1, 3.0
    )
    ecog = [a for a in alerts if a["alert_type"] == "ECOG_DECLINE"][0]
    assert ecog["severity"] == "HIGH"
    assert "terminal" not in ecog["message"].lower()


def test_child_pugh_c_alert_is_assessment_not_stage_assignment():
    alerts = BCLCAlertEngine().evaluate(
        "P1", "A", True, 0, "C", 1, 3.0
    )
    message = [a for a in alerts if a["alert_type"] == "SEVERE_LIVER_DYSFUNCTION"][0]["message"]
    assert "does not define BCLC-D" in message


def test_longitudinal_uses_core_staging_logic():
    assessment = HCCAssessment(date="2026-01-01", tumor_size_cm=7.0, tumor_count=1)
    assert bclc_stage(assessment) == "A"


def test_optn_downstaging_screen_current_boundaries():
    single = HCCAssessment(date="2026-01-01", tumor_size_cm=6.0, tumor_count=1, afp_ng_ml=100)
    assert optn_downstaging_screen(single)["screen_positive"] is True

    two = HCCAssessment(
        date="2026-01-01",
        tumor_size_cm=4.0,
        tumor_count=2,
        total_tumor_diameter_cm=7.0,
        afp_ng_ml=100,
    )
    assert optn_downstaging_screen(two)["screen_positive"] is True

    four_boundary = HCCAssessment(
        date="2026-01-01",
        tumor_size_cm=3.0,
        tumor_count=4,
        total_tumor_diameter_cm=8.0,
        afp_ng_ml=100,
    )
    assert optn_downstaging_screen(four_boundary)["screen_positive"] is False

    high_afp = HCCAssessment(date="2026-01-01", tumor_size_cm=6.0, tumor_count=1, afp_ng_ml=1001)
    assert optn_downstaging_screen(high_afp)["screen_positive"] is False
    assert "500" in optn_downstaging_screen(high_afp)["policy_note"]


def test_stratification_does_not_claim_transplant_eligibility():
    result = stratify_treatment("A", True, "A", 0, 3.0, 1)
    assert result["transplant_eligible"] is None
    assert result["transplant_evaluation_required"] is True


def test_dependency_free_simulator_runs():
    result = run_simulation(25, seed=1)
    assert sum(result["stage_counts"].values()) == 25
    assert result["invalid_input_checks_blocked"] == 1
