from src.analysis import get_cell_frequency_data
from src.statistics import compare_responders


def main() -> None:
    print("=== Part 2: Data Overview (First 5 rows) ===")
    df = get_cell_frequency_data()
    part2_df = df.rename(columns={"sample_id": "sample", "cell_type": "population"})
    display_cols = ["sample", "total_count", "population", "count", "percentage"]
    print(part2_df[display_cols].head().to_string(index=False))
    print(f"\nTotal rows processed: {len(df)}")

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


if __name__ == "__main__":
    main()
