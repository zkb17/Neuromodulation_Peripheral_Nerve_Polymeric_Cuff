"""
Statistical analysis script for:
  1) regressionData.csv (selectivity ex vivo and simulation)
  2) ElectrodeAgreement_SelectivityMagnitude_Scores_Top10.csv (electrode agreement)

Computes the statistics used in figure captions:

Figure Y — Ex vivo only, across nerves (max-per-config averaged within nerve):
  - For each fascicle:
      For each nerve: take max selectivity within each config (5 values), then average across configs.
    Then compare nerves:
      One-way ANOVA across nerves (p-value)
      Kruskal–Wallis across nerves (p-value)

Figure Z — Ex vivo vs Simulation, across waveform configurations:
  - For each fascicle:
      For each nerve×config: take max selectivity within that config for ExVivo and Simulation.
    Then compare modalities:
      Paired t-test across matched nerve×config pairs (p-value)

Figure W — Ex vivo vs Simulation, across nerves:
  - For each fascicle:
      For each nerve: max per config then average across configs (one value per nerve) for ExVivo and Simulation.
    Then compare modalities:
      Paired t-test across nerves (p-value)

Figure EA — Electrode agreement, across nerves:
  - For each fascicle:
      For each nerve: use all configuration-level electrode agreement values (configs ignored as factors).
    Then compare nerves:
      One-way ANOVA across nerves (p-value)
      Kruskal–Wallis across nerves (p-value)

Notes:
  - The “swap nerve 1 and 2 for display” does not change p-values.
  - Requires: pandas, numpy, scipy
"""

import numpy as np
import pandas as pd
from scipy import stats

# -------------------------
# User settings
# -------------------------
REGRESSION_CSV = "regressionData.csv"
ELECTRODE_AGREEMENT_CSV = "ElectrodeAgreement_SelectivityMagnitude_Scores_Top10.csv"

FASCICLES = ["Tibial", "Peroneal", "Sural"]
CONFIGS = [300, 303, 307, 305, 306]

# Optional: plotting/display order (stats invariant to order)
NERVE_ORDER_FOR_DISPLAY = [2, 1, 3, 4, 5]
NERVE_LABELS_FOR_DISPLAY = ["1", "2", "3", "4", "5"]

# regressionData.csv columns
COL_NERVE = "AnimalID"
COL_FASC = "Fascicle"
COL_CFG = "configurationID"
COL_EX = "Fascicle_Selectivity_ExVivo"
COL_SIM = "Fascicle_Selectivity_Simulation"

# electrode agreement columns
EA_COL_NERVE = "AnimalID"
EA_COL_FASC = "Fascicle"
EA_COL_SCORE = "ElectrodeAgreement"
EA_COL_CFG = "ConfigID"  # present but not used as a factor here


# -------------------------
# Helpers
# -------------------------
def load_regression(csv_file: str) -> pd.DataFrame:
    df = pd.read_csv(csv_file)
    df = df.dropna(subset=[COL_NERVE, COL_FASC, COL_CFG, COL_EX, COL_SIM]).copy()
    df = df[df[COL_FASC].isin(FASCICLES)]
    df = df[df[COL_CFG].isin(CONFIGS)]
    return df


def load_electrode_agreement(csv_file: str) -> pd.DataFrame:
    df = pd.read_csv(csv_file)
    df = df.dropna(subset=[EA_COL_NERVE, EA_COL_FASC, EA_COL_SCORE]).copy()
    df = df[df[EA_COL_FASC].isin(FASCICLES)]
    return df


def max_per_config(df: pd.DataFrame, value_col: str) -> pd.DataFrame:
    """One row per (nerve, fascicle, config): Value = max(value_col) within that nerve×fasc×config."""
    return (
        df.groupby([COL_NERVE, COL_FASC, COL_CFG])[value_col]
          .max()
          .reset_index()
          .rename(columns={value_col: "Value"})
    )


def average_across_configs(max_table: pd.DataFrame) -> pd.DataFrame:
    """One row per (nerve, fascicle): MeanAcrossConfigs = mean(Value across configs)."""
    return (
        max_table.groupby([COL_NERVE, COL_FASC])["Value"]
                 .mean()
                 .reset_index()
                 .rename(columns={"Value": "MeanAcrossConfigs"})
    )


def paired_ttest(x: np.ndarray, y: np.ndarray):
    mask = ~np.isnan(x) & ~np.isnan(y)
    if mask.sum() < 2:
        return np.nan, np.nan
    t, p = stats.ttest_rel(x[mask], y[mask])
    return float(t), float(p)


def one_way_anova(groups: list[np.ndarray]):
    groups = [g[~np.isnan(g)] for g in groups if np.sum(~np.isnan(g)) > 0]
    if len(groups) < 2:
        return np.nan, np.nan
    f, p = stats.f_oneway(*groups)
    return float(f), float(p)


def kruskal(groups: list[np.ndarray]):
    groups = [g[~np.isnan(g)] for g in groups if np.sum(~np.isnan(g)) > 0]
    if len(groups) < 2:
        return np.nan, np.nan
    h, p = stats.kruskal(*groups)
    return float(h), float(p)


def header(title: str):
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)


# -------------------------
# Figure Y — Ex vivo only, across nerves (max-per-config averaged)
# -------------------------
def stats_exvivo_across_nerves_max_per_config(df: pd.DataFrame):
    max_ex = max_per_config(df, COL_EX)
    mean_ex = average_across_configs(max_ex)

    nerves = sorted(mean_ex[COL_NERVE].unique())

    header("Figure Y — Ex vivo only, across nerves (max-per-config averaged within nerve)")

    for fasc in FASCICLES:
        sub = mean_ex[mean_ex[COL_FASC] == fasc]
        groups = [sub.loc[sub[COL_NERVE] == n, "MeanAcrossConfigs"].to_numpy(dtype=float) for n in nerves]

        F, pA = one_way_anova(groups)
        H, pK = kruskal(groups)

        print(f"\n{fasc}:")
        print(f"  One-way ANOVA p = {pA:.6g}" + (f"  (F = {F:.4g})" if not np.isnan(F) else ""))
        print(f"  Kruskal–Wallis p = {pK:.6g}" + (f"  (H = {H:.4g})" if not np.isnan(H) else ""))


# -------------------------
# Figure Z — Ex vivo vs simulation, across configs (paired across nerve×config)
# -------------------------
def stats_exvivo_vs_sim_across_configs(df: pd.DataFrame):
    max_ex = max_per_config(df, COL_EX).rename(columns={"Value": "ExVivo"})
    max_sim = max_per_config(df, COL_SIM).rename(columns={"Value": "Sim"})

    merged = pd.merge(max_ex, max_sim, on=[COL_NERVE, COL_FASC, COL_CFG], how="inner")

    header("Figure Z — Ex vivo vs Simulation, across configs (paired across nerve×config pairs)")

    for fasc in FASCICLES:
        sub = merged[merged[COL_FASC] == fasc]
        x = sub["ExVivo"].to_numpy(dtype=float)
        y = sub["Sim"].to_numpy(dtype=float)

        t, p = paired_ttest(x, y)
        n_pairs = int(np.sum(~np.isnan(x) & ~np.isnan(y)))
        print(f"\n{fasc}: paired t-test p = {p:.6g}" + (f"  (n = {n_pairs})" if not np.isnan(p) else ""))


# -------------------------
# Figure W — Ex vivo vs simulation, across nerves (max-per-config averaged; paired across nerves)
# -------------------------
def stats_exvivo_vs_sim_across_nerves_max_per_config(df: pd.DataFrame):
    max_ex = max_per_config(df, COL_EX)
    max_sim = max_per_config(df, COL_SIM)

    mean_ex = average_across_configs(max_ex).rename(columns={"MeanAcrossConfigs": "ExVivo"})
    mean_sim = average_across_configs(max_sim).rename(columns={"MeanAcrossConfigs": "Sim"})

    merged = pd.merge(mean_ex, mean_sim, on=[COL_NERVE, COL_FASC], how="inner").sort_values(COL_NERVE)

    header("Figure W — Ex vivo vs Simulation, across nerves (max-per-config averaged; paired across nerves)")

    for fasc in FASCICLES:
        sub = merged[merged[COL_FASC] == fasc]
        x = sub["ExVivo"].to_numpy(dtype=float)
        y = sub["Sim"].to_numpy(dtype=float)

        t, p = paired_ttest(x, y)
        n = int(np.sum(~np.isnan(x) & ~np.isnan(y)))
        print(f"\n{fasc}: paired t-test p = {p:.6g}" + (f"  (n = {n})" if not np.isnan(p) else ""))


# -------------------------
# Figure EA — Electrode agreement across nerves
# -------------------------
def stats_electrode_agreement_across_nerves(df_ea: pd.DataFrame):
    nerves = sorted(df_ea[EA_COL_NERVE].unique())

    header("Figure EA — Electrode agreement, across nerves (configs not treated as a factor)")

    for fasc in FASCICLES:
        sub = df_ea[df_ea[EA_COL_FASC] == fasc]
        groups = [sub.loc[sub[EA_COL_NERVE] == n, EA_COL_SCORE].to_numpy(dtype=float) for n in nerves]

        F, pA = one_way_anova(groups)
        H, pK = kruskal(groups)

        print(f"\n{fasc}:")
        print(f"  One-way ANOVA p = {pA:.6g}" + (f"  (F = {F:.4g})" if not np.isnan(F) else ""))
        print(f"  Kruskal–Wallis p = {pK:.6g}" + (f"  (H = {H:.4g})" if not np.isnan(H) else ""))


# -------------------------
# Main
# -------------------------
if __name__ == "__main__":
    reg = load_regression(REGRESSION_CSV)
    ea = load_electrode_agreement(ELECTRODE_AGREEMENT_CSV)

    # Selectivity stats
    stats_exvivo_across_nerves_max_per_config(reg)     # Figure Y
    stats_exvivo_vs_sim_across_configs(reg)            # Figure Z
    stats_exvivo_vs_sim_across_nerves_max_per_config(reg)  # Figure W

    # Electrode agreement stats
    stats_electrode_agreement_across_nerves(ea)        # Figure EA

    print("\nDone.")
