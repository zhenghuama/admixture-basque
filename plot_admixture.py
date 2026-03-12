import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re


def plot_admixture_results(
    fig_title,
    sampleinfo_path,
    fam_path,
    sample_name,
    pop_code_name,
    q_path_template,
    k_values,
    figsize=(10, 8),
):
    sampleinfo = pd.read_csv(sampleinfo_path, sep="\t")
    sample_to_pop = dict(zip(sampleinfo[sample_name], sampleinfo[pop_code_name]))

    with open(fam_path, "r") as f:
        samples = [line.split()[0] for line in f.readlines()]
    pops = [sample_to_pop.get(item, "NA") for item in samples]

    fig = plt.figure(figsize=figsize)

    for plotind, K in enumerate(k_values, start=1):
        ax = fig.add_subplot(3, 1, plotind)

        q_path = q_path_template.format(K=K)
        data = pd.read_csv(q_path, sep=r"\s+", header=None)

        cols = list(data.columns)
        data["sample"] = samples
        data["pop"] = pops
        data = data.sort_values(["pop"] + cols)

        data[cols].plot.bar(stacked=True, ax=ax, width=1, legend=False)

        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.yaxis.set_ticks_position("left")
        ax.xaxis.set_ticks_position("bottom")
        ax.set_title(f"K = {K}")

        # Only label first sample in each population
        xticklabels = []
        currpop = ""
        for i in range(data.shape[0]):
            pop_i = data["pop"].values[i]
            if pop_i == currpop:
                xticklabels.append("")
            else:
                xticklabels.append(pop_i)
                currpop = pop_i

        ax.set_xticks(range(len(xticklabels)))
        ax.set_xticklabels(xticklabels, rotation=90)

    fig.suptitle(fig_title)
    fig.tight_layout()
    return fig

# This Helper function is written by LLM
# It reads the reference result file and extracts the log-likelihood value from the line that starts with "Loglikelihood:"
def _extract_reference_loglikelihood(reference_result_path):
    pattern = re.compile(r"^Loglikelihood:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*$")

    with open(reference_result_path, "r") as f:
        for line in f:
            match = pattern.match(line.strip())
            if match is not None:
                return float(match.group(1))

    raise ValueError(
        f"Could not find a summary 'Loglikelihood:' line in {reference_result_path}."
    )


def compare_em_squarem(
    fig_title = "EM vs SQUAREM Convergence",
    reference_result_path = "data_and_results/step1/reference_results/step1_admixture.3.0",
    em_npy_path = "data_and_results/step1/results_em/step1_admixture.pruned.3.npy",
    squarem_npy_path = "data_and_results/step1/results_squarem/step1_admixture.pruned.3.npy",
    figsize=(10, 5),
):
    em_llh = np.load(em_npy_path)
    squarem_llh = np.load(squarem_npy_path)
    reference_llh = _extract_reference_loglikelihood(reference_result_path)

    em_iters = np.arange(1, em_llh.size + 1)
    squarem_iters = np.arange(1, squarem_llh.size + 1)

    fig, ax_llh = plt.subplots(figsize=figsize)

    ax_llh.plot(em_iters, em_llh, label=f"EM ({em_llh.size} iters)", linewidth=2)
    ax_llh.plot(
        squarem_iters,
        squarem_llh,
        label=f"SQUAREM ({squarem_llh.size} iters)",
        linewidth=2,
    )
    ax_llh.axhline(
        reference_llh,
        color="black",
        linestyle=":",
        linewidth=2,
        label="ADMIXTURE Convergence result",
    )
    ax_llh.set_title(fig_title)
    ax_llh.set_xlabel("Iteration")
    ax_llh.set_ylabel("Log-likelihood")
    ax_llh.grid(True, alpha=0.3)
    ax_llh.legend()

    fig.tight_layout()
    return fig
