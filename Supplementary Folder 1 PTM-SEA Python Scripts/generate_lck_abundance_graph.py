"""
====================================================================================================
Kenneth P. Callahan


7 September 2021

====================================================================================================
Python >= 3.8.10

generate_lck_abundance_graph.py
====================================================================================================
This python file accompanies the manuscript 

CD19-CAR T cells inhibit B-cell activating phosphorylation in target cells
by Griffith et al.

and generates Supplementary Figure 11.
====================================================================================================
"""

####################################################################################################
#
#   Importables

import pandas as pd
import matplotlib.pyplot as plt
from py_scripts.helpers.homebrew_stats import mean, standard_deviation, pairwise_t

#
#
####################################################################################################
#
#   Functions

def add_errorbar(mpl_axes, x_pos, y_pos, std, color = "grey", x_offset = 0.05, transparency = 0.75):
    """
    Adds error bars to MatPlotLib axes object. Modifies the mpl_axes input, returns None.

    mpl_axes -> matplotlib axes object
    x_pos    -> position of the center of the error bars on the x axis
    y_pos    -> position of the center of the error bars on the y axis
    std      -> standard deviation of the data (or vertical offset for error bars)
    color    -> The color of the error bars (Default: "grey")
    x_offset -> Size of the middle of the error bars (Default: 0.05)
    """
    # Plot the vertical middle bar
    mpl_axes.plot([x_pos, x_pos], [y_pos + std, y_pos - std], color = color, alpha = transparency)
    # Plot the horizontal middle bar
    mpl_axes.plot([x_pos - x_offset, x_pos + x_offset], [y_pos, y_pos], color = color, alpha = transparency)
    # Plot the horizontal top bar
    mpl_axes.plot([x_pos - 0.75*x_offset, x_pos + 0.75*x_offset], [y_pos + std, y_pos + std], color = color, alpha = transparency)
    # Plot the horizontal bottom bar
    mpl_axes.plot([x_pos - 0.75*x_offset, x_pos + 0.75*x_offset], [y_pos - std, y_pos - std], color = color, alpha = transparency)
    return None

def make_sigstrings(pvals, symbol = "*"):
    """
    """
    sigstrings = []
    for p in pvals:
        if p >= 0.05:
            sigstrings.append(f"n.s.")
        elif 0.01 <= p < 0.05:
            sigstrings.append(f"{symbol}")
        elif 0.001 <= p < 0.01:
            sigstrings.append(f"{symbol}{symbol}")
        else:
            sigstrings.append(f"{symbol}{symbol}{symbol}")
    return sigstrings

#
#
####################################################################################################
#
#  main() function

def main():
    """
    """
    file = pd.read_excel("Supplementary Table 6 Lck Western Quantification.xls")
    signal = list(file["Unnamed: 3"].to_numpy())
    tps = signal[1:17]
    lck = signal[20:]
    normalized = [lck[i] / tps[i] for i in range(len(tps))]
    samples = [normalized[i*4:(i+1)*4] for i in range(4)]
    # The last group of data is just noise, so we omit it.
    del samples[-1]
    samples = [("JE6",samples[0]), ("CD19-CAR T cells",samples[1]), ("CD19$^{HI}$ Raji B cells",samples[2])]
    lck_ttest = pairwise_t(*samples,
                            t_type = "s")
    print(lck_ttest)
    lck_ttest.to_excel("lck_ttest_results.xlsx")
    sigs = make_sigstrings(list(lck_ttest['pvalue'].astype(float).to_numpy()))
    xs = [[i-0.2+0.05*j for j in range(4)] for i in range(1,4)]
    centers = [i-0.2 + (0.05*3/2) for i in range(1,4)]
    maxs = [max(samples[i][1]) for i in range(len(samples))]
    means = [mean(item[1]) for item in samples]
    stds = [standard_deviation(item[1]) for item in samples]
    mean_up = [means[i] + stds[i] for i in range(len(means))]
    mean_down = [means[i] - stds[i] for i in range(len(means))]
    fig, ax = plt.subplots()
    for i in range(3):
        ax.scatter(xs[i], samples[i][1], edgecolors = "black")
        add_errorbar(ax, centers[i], means[i], stds[i])
    ax.set_ylim(0,.23)
    ax.plot([centers[0]+0.005, centers[0]+0.005], [maxs[0]+0.018, maxs[0]+0.03], color = "black")
    ax.plot([centers[1]-0.005, centers[1]-0.005], [maxs[1]+0.018, maxs[0]+0.03], color = "black")
    ax.plot([centers[0]+0.005, centers[1]-0.005], [maxs[0]+0.03,maxs[0]+0.03], color = "black")
    ax.plot([centers[0]-0.005, centers[0]-0.005], [maxs[0]+0.018, maxs[0]+0.05], color = "black")
    ax.plot([centers[2]+0.005, centers[2]+0.005], [maxs[2] + 0.018, maxs[0] + 0.05], color = "black")
    ax.plot([centers[0]-0.005, centers[2]+0.005], [maxs[0] + 0.05, maxs[0] + 0.05], color = "black")
    ax.plot([centers[1]+0.005, centers[1]+0.005], [maxs[1] + 0.018, maxs[1] + 0.03], color = "black")
    ax.plot([centers[2]-0.005, centers[2] - 0.005], [maxs[2] +0.018, maxs[1] + 0.03], color = "black")
    ax.plot([centers[1] + 0.005, centers[2] - 0.005], [maxs[1] + 0.03, maxs[1] + 0.03], color = "black")
    ax.text((centers[0] + centers[1])/2, maxs[0] + 0.032, sigs[0])
    ax.text((centers[0] + centers[2])/2, maxs[0] + 0.052, sigs[1])
    ax.text((centers[1] + centers[2])/2, maxs[1] + 0.032, sigs[2])
    ax.set_title("Lck Abundance")
    ax.set_ylabel("Lck signal / Total Protein signal")
    ax.set_xticks(centers)
    ax.set_xticklabels([samples[0][0], samples[1][0], samples[2][0]])#, rotation = 0, ha = "center", rotation_mode = "anchor")
    plt.savefig("lck_abundance.pdf", bbox_inches = "tight")
    plt.show()
    return None

#
#
#################################################################################################

    
if __name__ == "__main__":
    main()