import json
from pathlib import Path
from typing import cast

from src.analysis import get_part2_frequency_table
from src.queries import get_subset_stats
from src.statistics import compare_responders


def main() -> None:
    output_dir = Path("outputs")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== Part 2: Data Overview (First 5 rows) ===")
    part2_df = get_part2_frequency_table()
    print(part2_df.head().to_string(index=False))
    print(f"\nTotal rows processed: {len(part2_df)}")
    part2_df.to_csv(output_dir / "part2_frequency_table.csv", index=False)

    print("\n=== Part 3: Statistical Analysis (Responders vs Non-Responders) ===")
    print("Condition: Melanoma, Treatment: Miraclib, Sample: PBMC")

    stats_df, _, summary = compare_responders(time_filter="baseline_only", unit="subject")
    print(
        f"Test: {summary['test_label']} | Correction: {summary['correction_label']} | "
        f"Unit: {summary['unit']} | Metric: {summary['metric']} | {summary['bootstrap_ci']}"
    )
    print(
        stats_df[
            [
                "cell_type",
                "n_yes",
                "n_no",
                "median_diff",
                "direction",
                "ci_95_low",
                "ci_95_high",
                "effect",
                "p_value",
                "q_value",
                "significant",
                "avg_responder",
                "avg_non_responder",
            ]
        ].to_string(index=False)
    )
    stats_df.to_csv(output_dir / "part3_stats.csv", index=False)

    print("\nInterpretation:")
    for row in stats_df.to_dict(orient="records"):
        is_significant = bool(row.get("significant", False))
        q_value_obj = row.get("q_value")
        q_value = float(q_value_obj) if isinstance(q_value_obj, (int, float)) else None
        if is_significant and q_value is not None:
            print(f"-> SIGNIFICANT difference found in {row['cell_type']} (q={q_value:.4f})")
            direction = str(row.get("direction", "undetermined"))
            ci_low = row.get("ci_95_low")
            ci_high = row.get("ci_95_high")
            print(f"   Direction: {direction}")
            if isinstance(ci_low, (int, float)) and isinstance(ci_high, (int, float)):
                print(f"   95% bootstrap CI: [{float(ci_low):.4f}, {float(ci_high):.4f}]")

    part4 = get_subset_stats(
        condition="melanoma",
        treatment="miraclib",
        sample_type="PBMC",
        time_filter="baseline_only",
    )
    n_projects = cast(int, part4["n_projects"])
    n_samples = cast(int, part4["n_samples"])
    n_subjects = cast(int, part4["n_subjects"])
    avg_b_raw = part4["avg_b_cell_male_responders"]
    part4_summary = {
        "n_projects": n_projects,
        "n_samples": n_samples,
        "n_subjects": n_subjects,
        "avg_b_cell_male_responders": (
            round(float(avg_b_raw), 2)
            if isinstance(avg_b_raw, float)
            else None
        ),
    }
    (output_dir / "part4_summary.json").write_text(
        json.dumps(part4_summary, indent=2),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
