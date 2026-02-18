import pandas as pd

from src.analysis import get_cell_frequency_data
from src.statistics import compare_responders


def main():
    print("=== Part 2: Data Overview (First 5 rows) ===")
    df = get_cell_frequency_data()
    display_cols = ["sample_id", "total_count", "cell_type", "count", "percentage"]
    print(df[display_cols].head().to_string(index=False))
    print(f"\nTotal rows processed: {len(df)}")

    print("\n=== Part 3: Statistical Analysis (Responders vs Non-Responders) ===")
    print("Condition: Melanoma, Treatment: Miraclib, Sample: PBMC")

    stats_df, _ = compare_responders()
    print(
        stats_df[
            ["cell_type", "p_value", "significant", "avg_responder", "avg_non_responder"]
        ].to_string(index=False)
    )

    print("\nInterpretation:")
    for row in stats_df.to_dict(orient="records"):
        is_significant = bool(row.get("significant", False))
        p_value_obj = row.get("p_value")
        p_value = float(p_value_obj) if isinstance(p_value_obj, (int, float)) else None
        if is_significant and p_value is not None:
            print(f"-> SIGNIFICANT difference found in {row['cell_type']} (p={p_value:.4f})")
            avg_responder = row.get("avg_responder")
            avg_non_responder = row.get("avg_non_responder")
            if isinstance(avg_responder, (int, float)) and isinstance(avg_non_responder, (int, float)) and float(avg_responder) > float(avg_non_responder):
                print("   (Higher in Responders)")
            else:
                print("   (Higher in Non-Responders)")


if __name__ == "__main__":
    main()
