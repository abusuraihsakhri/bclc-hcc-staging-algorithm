#!/usr/bin/env python3
"""Tests for BCLC HCC Staging Algorithm - 20 real clinical tests."""

import json
import pytest
from bclc_staging import (
    stage_bclc,
    process_batch,
    _check_milan,
    _parse_bool,
    BCLC_STAGES,
)


# ---------------------------------------------------------------------------
# Stage 0 - Very Early
# ---------------------------------------------------------------------------

class TestStage0:
    def test_stage_0_single_small(self):
        """Single ≤2cm, Child-Pugh A, ECOG 0 → Stage 0."""
        res = stage_bclc(tumor_count=1, tumor_size_cm=1.5, child_pugh_class="A", ecog_ps=0)
        assert res["bclc_stage"] == "0"
        assert res["stage_name"] == "Very Early"

    def test_stage_0_boundary_2cm(self):
        """Single exactly 2cm → Stage 0."""
        res = stage_bclc(tumor_count=1, tumor_size_cm=2.0, child_pugh_class="A", ecog_ps=0)
        assert res["bclc_stage"] == "0"

    def test_stage_0_treatment(self):
        res = stage_bclc(tumor_count=1, tumor_size_cm=1.5, child_pugh_class="A", ecog_ps=0)
        assert "Ablation" in res["treatment_allocation"] or "Resection" in res["treatment_allocation"]


# ---------------------------------------------------------------------------
# Stage A - Early
# ---------------------------------------------------------------------------

class TestStageA:
    def test_stage_a_single_3cm(self):
        """Single 3cm tumor, Child-Pugh A, ECOG 0 → Stage A."""
        res = stage_bclc(tumor_count=1, tumor_size_cm=3.0, child_pugh_class="A", ecog_ps=0)
        assert res["bclc_stage"] == "A"
        assert res["stage_name"] == "Early"

    def test_stage_a_multi_small(self):
        """3 nodules each ≤3cm → Stage A."""
        res = stage_bclc(tumor_count=3, tumor_size_cm=2.5, child_pugh_class="A", ecog_ps=0)
        assert res["bclc_stage"] == "A"

    def test_stage_a_single_up_to_5cm(self):
        """Single tumor up to 5cm, Child-Pugh A → Stage A."""
        res = stage_bclc(tumor_count=1, tumor_size_cm=4.5, child_pugh_class="A", ecog_ps=0)
        assert res["bclc_stage"] == "A"

    def test_stage_a_treatment(self):
        res = stage_bclc(tumor_count=1, tumor_size_cm=3.0, child_pugh_class="A", ecog_ps=0)
        assert "Resection" in res["treatment_allocation"] or "Transplant" in res["treatment_allocation"]


# ---------------------------------------------------------------------------
# Stage B - Intermediate
# ---------------------------------------------------------------------------

class TestStageB:
    def test_stage_b_multinodular(self):
        """Multinodular (>3), Child-Pugh A, ECOG 0 → Stage B."""
        res = stage_bclc(tumor_count=5, tumor_size_cm=3.0, child_pugh_class="A", ecog_ps=0)
        assert res["bclc_stage"] == "B"
        assert res["stage_name"] == "Intermediate"

    def test_stage_b_child_pugh_b(self):
        """Child-Pugh B with small tumor → Stage B."""
        res = stage_bclc(tumor_count=1, tumor_size_cm=2.0, child_pugh_class="B", ecog_ps=0)
        assert res["bclc_stage"] == "B"

    def test_stage_b_treatment(self):
        res = stage_bclc(tumor_count=5, tumor_size_cm=3.0, child_pugh_class="A", ecog_ps=0)
        assert "TACE" in res["treatment_allocation"]


# ---------------------------------------------------------------------------
# Stage C - Advanced
# ---------------------------------------------------------------------------

class TestStageC:
    def test_stage_c_portal_invasion(self):
        """Portal vein invasion → Stage C."""
        res = stage_bclc(
            tumor_count=1, tumor_size_cm=3.0, child_pugh_class="A", ecog_ps=0,
            portal_vein_invasion=True,
        )
        assert res["bclc_stage"] == "C"
        assert res["stage_name"] == "Advanced"

    def test_stage_c_extrahepatic(self):
        """Extrahepatic spread → Stage C."""
        res = stage_bclc(
            tumor_count=1, tumor_size_cm=3.0, child_pugh_class="A", ecog_ps=0,
            extrahepatic_spread=True,
        )
        assert res["bclc_stage"] == "C"

    def test_stage_c_lymph_nodes(self):
        """Lymph node metastasis → Stage C."""
        res = stage_bclc(
            tumor_count=1, tumor_size_cm=3.0, child_pugh_class="A", ecog_ps=0,
            lymph_node_metastasis=True,
        )
        assert res["bclc_stage"] == "C"

    def test_stage_c_ecog_1(self):
        """ECOG 1 with Child-Pugh A → Stage C."""
        res = stage_bclc(tumor_count=1, tumor_size_cm=3.0, child_pugh_class="A", ecog_ps=1)
        assert res["bclc_stage"] == "C"

    def test_stage_c_treatment(self):
        res = stage_bclc(
            tumor_count=1, tumor_size_cm=3.0, child_pugh_class="A", ecog_ps=1,
        )
        assert "Systemic" in res["treatment_allocation"]


# ---------------------------------------------------------------------------
# Stage D - Terminal
# ---------------------------------------------------------------------------

class TestStageD:
    def test_stage_d_child_pugh_c(self):
        """Child-Pugh C → Stage D regardless of tumor."""
        res = stage_bclc(tumor_count=1, tumor_size_cm=1.0, child_pugh_class="C", ecog_ps=0)
        assert res["bclc_stage"] == "D"
        assert res["stage_name"] == "Terminal"

    def test_stage_d_ecog_3(self):
        """ECOG 3 → Stage D."""
        res = stage_bclc(tumor_count=1, tumor_size_cm=1.0, child_pugh_class="A", ecog_ps=3)
        assert res["bclc_stage"] == "D"

    def test_stage_d_treatment(self):
        res = stage_bclc(tumor_count=1, tumor_size_cm=1.0, child_pugh_class="C", ecog_ps=0)
        assert "Supportive" in res["treatment_allocation"]


# ---------------------------------------------------------------------------
# Milan criteria
# ---------------------------------------------------------------------------

class TestMilanCriteria:
    def test_milan_single_small(self):
        assert _check_milan(1, 3.0) is True

    def test_milan_single_5cm(self):
        assert _check_milan(1, 5.0) is True

    def test_milan_single_over_5cm(self):
        assert _check_milan(1, 6.0) is False

    def test_milan_multi_3_nodules_3cm(self):
        assert _check_milan(3, 3.0) is True

    def test_milan_multi_4_nodules(self):
        assert _check_milan(4, 2.0) is False

    def test_milan_multi_over_3cm(self):
        assert _check_milan(2, 4.0) is False


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

class TestValidation:
    def test_invalid_child_pugh(self):
        with pytest.raises(ValueError):
            stage_bclc(child_pugh_class="D")

    def test_invalid_ecog(self):
        with pytest.raises(ValueError):
            stage_bclc(ecog_ps=5)

    def test_result_has_all_fields(self):
        res = stage_bclc(tumor_count=1, tumor_size_cm=2.0, child_pugh_class="A", ecog_ps=0)
        assert "bclc_stage" in res
        assert "treatment_allocation" in res
        assert "milan_criteria_eligible" in res
        assert "five_year_survival_pct" in res
        assert "clinical_recommendation" in res
        assert "classification" in res


# ---------------------------------------------------------------------------
# Batch processing
# ---------------------------------------------------------------------------

class TestBatch:
    def test_batch_basic(self, tmp_path):
        csv_in = tmp_path / "in.csv"
        csv_out = tmp_path / "out.csv"
        csv_in.write_text(
            "tumor_count,tumor_size_cm,child_pugh_class,ecog_ps\n"
            "1,1.5,A,0\n"
            "5,3.0,B,0\n",
            encoding="utf-8",
        )
        count = process_batch(str(csv_in), str(csv_out))
        assert count == 2
        assert csv_out.exists()
        content = csv_out.read_text(encoding="utf-8")
        assert "bclc_stage" in content
        assert "treatment_allocation" in content


# ---------------------------------------------------------------------------
# Parse bool helper
# ---------------------------------------------------------------------------

class TestParseBool:
    def test_parse_bool_true_values(self):
        assert _parse_bool("true") is True
        assert _parse_bool("1") is True
        assert _parse_bool("yes") is True

    def test_parse_bool_false_values(self):
        assert _parse_bool("false") is False
        assert _parse_bool("0") is False
        assert _parse_bool("no") is False
