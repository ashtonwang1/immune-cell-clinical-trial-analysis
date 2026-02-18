from typing import cast

import pandas as pd
from scipy import stats

from src.analysis import apply_clr_transform, get_filtered_data, prepare_unit_level_data


def _bh_fdr_adjust(p_values: list[float | None]) -> list[float | None]:
    indexed = [(idx, p) for idx, p in enumerate(p_values) if p is not None]
    if not indexed:
        return [None] * len(p_values)

    indexed.sort(key=lambda x: x[1])
    m = len(indexed)
    adjusted = [0.0] * m

    min_so_far = 1.0
    for rank in range(m - 1, -1, -1):
        _, p = indexed[rank]
        q = (p * m) / (rank + 1)
        min_so_far = min(min_so_far, q)
        adjusted[rank] = min(min_so_far, 1.0)

    output: list[float | None] = [None] * len(p_values)
    for rank, (orig_idx, _) in enumerate(indexed):
        output[orig_idx] = float(adjusted[rank])
    return output


def _cliffs_delta(group_yes: list[float], group_no: list[float]) -> float:
    total = len(group_yes) * len(group_no)
    if total == 0:
        return 0.0
    gt = 0
    lt = 0
    for a in group_yes:
        for b in group_no:
            if a > b:
                gt += 1
            elif a < b:
                lt += 1
    return (gt - lt) / total


def compare_responders(
    condition: str = "melanoma",
    treatment: str = "miraclib",
    sample_type: str = "PBMC",
    time_filter: str = "all",
    unit: str = "sample",
    metric: str = "percentage",
    transform: str = "none",
    test: str = "mannwhitney",
    correction: str = "bh_fdr",
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, str]]:
    df = get_filtered_data(condition, treatment, sample_type, time_filter=time_filter)
    plot_df = prepare_unit_level_data(df, unit=unit, metric=metric)
    if transform == "clr":
        plot_df = apply_clr_transform(plot_df)

    results = []
    cell_types = sorted(set(plot_df["cell_type"].dropna().astype(str).tolist()))

    for cell in cell_types:
        subset = plot_df[plot_df["cell_type"] == cell]
        group_yes = subset.loc[subset["response"] == "yes", "metric_value"].tolist()
        group_no = subset.loc[subset["response"] == "no", "metric_value"].tolist()

        n_yes = len(group_yes)
        n_no = len(group_no)

        if group_yes and group_no:
            if test == "welch_t":
                test_result = stats.ttest_ind(group_yes, group_no, equal_var=False)
                p_value = cast(float, test_result[1])
                stat_score = cast(float, test_result[0])
                effect = float(pd.Series(group_yes).mean() - pd.Series(group_no).mean())
                effect_label = "mean_diff"
            else:
                stat, p_val = stats.mannwhitneyu(group_yes, group_no, alternative="two-sided")
                p_value = float(p_val)
                stat_score = float(stat)
                effect = float((2 * stat_score) / (n_yes * n_no) - 1)
                effect_label = "rank_biserial"
            cliffs = _cliffs_delta(group_yes, group_no)
        else:
            p_value = None
            stat_score = None
            effect = None
            cliffs = None
            effect_label = "rank_biserial" if test == "mannwhitney" else "mean_diff"

        results.append(
            {
                "cell_type": cell,
                "n_yes": n_yes,
                "n_no": n_no,
                "p_value": p_value,
                "stat_score": stat_score,
                "median_yes": float(pd.Series(group_yes).median()) if group_yes else None,
                "median_no": float(pd.Series(group_no).median()) if group_no else None,
                "effect": effect,
                "effect_label": effect_label,
                "cliffs_delta": cliffs,
                "avg_responder": float(pd.Series(group_yes).mean()) if group_yes else None,
                "avg_non_responder": float(pd.Series(group_no).mean()) if group_no else None,
            }
        )

    stats_df = pd.DataFrame(results)
    if len(stats_df) == 0:
        summary = {
            "test_label": "Welch t-test" if test == "welch_t" else "Mann-Whitney U",
            "correction_label": "None" if correction == "none" else "BH-FDR",
            "unit": unit,
            "metric": metric,
            "transform_label": "CLR" if transform == "clr" else "Raw",
        }
        return stats_df, plot_df, summary

    if correction == "none":
        stats_df["q_value"] = stats_df["p_value"]
    else:
        stats_df["q_value"] = _bh_fdr_adjust(stats_df["p_value"].tolist())

    stats_df["significant"] = stats_df["q_value"].apply(lambda q: bool(q is not None and q < 0.05))
    stats_df = stats_df.sort_values(by=["q_value", "p_value"], na_position="last").reset_index(drop=True)

    summary = {
        "test_label": "Welch t-test" if test == "welch_t" else "Mann-Whitney U",
        "correction_label": "None" if correction == "none" else "BH-FDR",
        "unit": unit,
        "metric": metric,
        "transform_label": "CLR" if transform == "clr" else "Raw",
    }

    return stats_df, plot_df, summary
