"""Read saved results for the Claude handoff; never refit or change source outputs."""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("evidence_snapshot.json"))
    args = parser.parse_args()
    root = args.project_root.resolve()
    manifest = {}

    def source(relative: str) -> Path:
        path = root / relative
        raw = path.read_bytes()
        manifest[relative] = {
            "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest(),
            "modified_utc": datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat(),
        }
        return path

    def csv(relative: str) -> pd.DataFrame:
        return pd.read_csv(source(relative))

    def metadata(relative: str) -> dict:
        return json.loads(source(relative).read_text(encoding="utf-8-sig"))

    base = "outputs/dynamic_social_unit_merge_gamm/daily_interaction_hurdle/"
    result = {
        "audit_utc": datetime.now(timezone.utc).isoformat(),
        "project_root": str(root),
        "scope": "Saved local metadata, CSVs and scripts only. No model reruns, source GPS validation, remote access or literature review.",
        "cohorts": {}, "meeting_input_audit": {},
    }
    for name, path in {
        "full_hourly_export": "outputs/membership_export_narrow/canonical_hourly_membership_narrow.metadata.json",
        "nightly_export": "outputs/membership_export_nightly/canonical_nightly_membership.metadata.json",
        "legacy_dynamic_status": "outputs/dynamic_social_unit_merge_gamm/proximity_status_dynamic_social_unit_metadata.json",
        "two_metre_merge_summary": "outputs/all_supported_2m_group_merges/all_supported_2m_group_merges_metadata.json",
        "copper_lilac": "outputs/copper_lilac_effort_corrected_integration/copper_lilac_effort_corrected_integration_metadata.json",
    }.items():
        result["cohorts"][name] = metadata(path)

    fit_columns = ["pair_key", "group_a", "group_b", "any_interaction",
                   "log1p_range_mean_centroid_distance_m_z", "log1p_prior_pair_interaction_days_z",
                   "log1p_group_size_total_z", "log1p_group_size_abs_diff_z", "dyad_daily_ndvi_mean_z"]
    for unit, filename in [("daily", "daily_interaction_hurdle_daily_event_rows.csv"),
                           ("weekly", "daily_interaction_hurdle_model_rows.csv")]:
        frame = csv(base + filename)
        fit = frame[fit_columns].replace([np.inf, -np.inf], np.nan).dropna()
        pair = frame[frame.pair_key.eq("Copper - Lilac")]
        result["meeting_input_audit"][unit] = {
            "candidate_rows": len(frame), "candidate_dyads": frame.pair_key.nunique(),
            "rows_with_zero_saved_eligible_bins": int(frame.eligible_2min_bins.eq(0).sum()),
            "copper_lilac_candidate_rows": len(pair),
            "copper_lilac_positive_rows": int(pair.any_interaction.sum()),
            "copper_lilac_saved_eligible_bins": int(pair.eligible_2min_bins.sum()),
            "ndvi_complete_fit_rows": len(fit), "ndvi_complete_fit_dyads": fit.pair_key.nunique(),
            "ndvi_complete_positive_rows": int(fit.any_interaction.sum()),
            "copper_lilac_ndvi_complete_fit_rows": int(fit.pair_key.eq("Copper - Lilac").sum()),
            "group_size_equals_sum_observed_animals_all_rows": bool(np.allclose(
                frame.group_size_total, frame.group_a_animals + frame.group_b_animals)),
        }

    for key, rel in {
        "meeting_model_summaries": base + "meeting_probability_summary.csv",
        "duration_model_summary": base + "event_duration_summary.csv",
        "integration_model_summary": base + "event_integration_summary.csv",
        "integration_metric_models": base + "integration_metric_summary.csv",
        "distance_window_comparison": base + "event_distance_window_comparison.csv",
        "weekly_gee_summary": "outputs/weekly_interaction_gamm/weekly_interaction_gamm_model_summary.csv",
        "copper_lilac_phase_ratios": "outputs/copper_lilac_effort_corrected_integration/copper_lilac_phase_integration_by_radius.csv",
        "copper_lilac_phase_dates": "outputs/copper_lilac_phase_alluvial_from_202403/copper_lilac_data_driven_phases.csv",
        "two_metre_dyad_summary": "outputs/all_supported_2m_group_merges/all_supported_2m_merge_dyad_summary.csv",
        "presplit_prediction": "outputs/presplit_cositting_prediction/prediction_summary.csv",
    }.items():
        result[key] = csv(rel).to_dict(orient="records")

    result["activity"] = {}
    for key, filename, column in [
        ("vedba", "event_activity_vedba_summary.csv", "event_vedba_mean"),
        ("gps_movement", "event_gps_movement_summary.csv", "event_gps_movement_mean_m_per_h"),
    ]:
        frame = csv(base + filename)
        result["activity"][key] = {
            "available_events": int(frame[column].notna().sum()),
            "pearson_r_with_integration_5m_fraction": frame[[column, "integration_5m_fraction"]].corr().iloc[0, 1],
        }
    old = csv("outputs/disperser_integration_calendar_time/disperser_event_outcomes_14d.csv")
    established = csv("outputs/established_dispersal_centrality/established_events.csv")
    cells = csv("outputs/established_dispersal_centrality/event_week_centrality_cells.csv")
    result["dispersal"] = {
        "earlier_events": len(old), "earlier_animals": old.animal_id.nunique(),
        "earlier_fates": old.fate.value_counts().to_dict(),
        "later_segments": len(established), "later_animals": established.animal_id.nunique(),
        "later_segment_labels": established.established_outcome.value_counts().to_dict(),
        "later_segments_with_centrality": cells.event_id.nunique(),
        "later_measured_by_label": cells.groupby("established_outcome").event_id.nunique().to_dict(),
    }
    for name in ["fit_daily_interaction_hurdle_model.py", "fit_weekly_interaction_gamm.py",
                 "analyze_copper_lilac_effort_corrected_integration.py", "analyze_presplit_cositting_predicts_split.py",
                 "compare_centrality_stayed_left.py", "analyze_established_dispersal_centrality.py",
                 "classify_disperser_event_outcomes.py", "plot_permanent_5m_colocation_cadence_filtered.py",
                 "docs/emerald_bronze_sparse_case_24AA11.md"]:
        source(name)

    result["source_manifest"] = manifest

    def clean(value):
        if isinstance(value, dict):
            return {str(k): clean(v) for k, v in value.items()}
        if isinstance(value, list):
            return [clean(v) for v in value]
        if isinstance(value, np.generic):
            return clean(value.item())
        if isinstance(value, float) and not np.isfinite(value):
            return None
        return value

    args.output.write_text(json.dumps(clean(result), indent=2, ensure_ascii=True, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Saved {args.output.resolve()}; inspected {len(manifest)} source files. Models were not rerun.")


if __name__ == "__main__":
    main()
