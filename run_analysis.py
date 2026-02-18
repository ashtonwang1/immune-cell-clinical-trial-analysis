import pandas as pd

from src.analysis import get_cell_frequency_data
from src.statistics import compare_responders


def main() -> None:
    print("=== Part 2: Data Overview (First 5 rows) ===")
    df = get_cell_frequency_data()
    display_cols = ["sample_id", "total_count", "cell_type", "count", "percentage"]
    print(df[display_cols].head().to_string(index=False))
    print(f"\nTotal rows processed: {len(df)}")

    print("\n=== Part 3: Statistical Analysis (Responders vs Non-Responders) ===")
    print("Condition: Melanoma, Treatment: Miraclib, Sample: PBMC")

    stats_df, _, summary = compare_responders()
    print(
        f"Test: {summary['test_label']} | Correction: {summary['correction_label']} | "
        f"Unit: {summary['unit']} | Metric: {summary['metric']}"
    )
    print(
        stats_df[
            [
                "cell_type",
                "n_yes",
                "n_no",
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
            avg_responder = row.get("avg_responder")
            avg_non_responder = row.get("avg_non_responder")
            if isinstance(avg_responder, (int, float)) and isinstance(avg_non_responder, (int, float)) and float(avg_responder) > float(avg_non_responder):
                print("   (Higher in Responders)")
            else:
                print("   (Higher in Non-Responders)")


if __name__ == "__main__":
    main()
