import pandas as pd
from scipy import stats

from src.analysis import get_filtered_data


def compare_responders(
    condition: str = "melanoma", treatment: str = "miraclib", sample_type: str = "PBMC"
) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = get_filtered_data(condition, treatment, sample_type)

    results = []
    cell_types = sorted(set(df["cell_type"].dropna().astype(str).tolist()))

    for cell in cell_types:
        subset = df[df["cell_type"] == cell]
        group_yes = subset.loc[subset["response"] == "yes", "percentage"].tolist()
        group_no = subset.loc[subset["response"] == "no", "percentage"].tolist()

        if group_yes and group_no:
            stat, p_val = stats.mannwhitneyu(group_yes, group_no, alternative="two-sided")
            p_value = float(p_val)
            stat_score = float(stat)
        else:
            p_value = None
            stat_score = None

        results.append(
            {
                "cell_type": cell,
                "p_value": p_value,
                "stat_score": stat_score,
                "significant": bool(p_value is not None and p_value < 0.05),
                "avg_responder": float(pd.Series(group_yes).mean()) if group_yes else None,
                "avg_non_responder": float(pd.Series(group_no).mean()) if group_no else None,
            }
        )

    return pd.DataFrame(results), df
