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
plt.rcParams["font.family"] = "sans-serif"
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
    correction_factor = [item / max(tps) for item in tps]
    normalized = [lck[i] / correction_factor[i] for i in range(len(tps))]
    samples = [normalized[i*4:(i+1)*4] for i in range(4)]
    # The last group of data is just noise, so we omit it.
    del samples[-1]
    samples = [("JE6",samples[0]), ("CD19-CAR T cells",samples[1]), ("Raji B cells",samples[2])]
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
    colours = ["blue", "cyan", "hotpink"]
    fig, ax = plt.subplots()
    for i in range(3):
        ax.scatter(xs[i], sorted(samples[i][1],reverse = True), color = colours[i], edgecolors = "black")
        add_errorbar(ax, centers[i], means[i], stds[i])
    #ax.set_ylim(0,.23)
    ax.plot([centers[1], centers[1]], [maxs[1]+180000, maxs[0]+400000], color = "cyan", linestyle = ":", alpha = 0.5)
    ax.plot([centers[0]+0.005, centers[1]-0.005], [maxs[0]+400000,maxs[0]+400000], color = "black", alpha = 0.5)
    ax.plot([centers[0], centers[0]], [maxs[0]+180000, maxs[0]+700000], color = "blue", linestyle = ":", alpha = 0.5)
    ax.plot([centers[2], centers[2]], [maxs[2] + 180000, maxs[0] + 700000], color = "hotpink", linestyle = ":", alpha = 0.5)
    ax.plot([centers[0], centers[2]], [maxs[0] + 700000, maxs[0] + 700000], color = "black", alpha = 0.5)
#    ax.plot([centers[1]+0.005, centers[1]+0.005], [maxs[1] + 700000, maxs[1] + 700000], color = "grey", linestyle = ":", alpha = 0.5)
#    ax.plot([centers[2]-0.005, centers[2] - 0.005], [maxs[2] +180000, maxs[1] + 400000], color = "grey", linestyle = ":", alpha = 0.5)
    ax.plot([centers[1], centers[2]], [maxs[1] + 400000, maxs[1] + 400000], color = "black", alpha = 0.5)
    ax.text((centers[0] + centers[1])/2, maxs[0] + 420000, sigs[0], font = "Arial", fontsize = 10, fontweight = "bold", ha = "center")
    ax.text((centers[0] + centers[2])/2, maxs[0] + 700000, sigs[1], font = "Arial", fontsize = 10, fontweight = "bold", ha = "center")
    ax.text((centers[1] + centers[2])/2, maxs[1] + 400000, sigs[2], font = "Arial", fontsize = 10, fontweight = "bold", ha = "center")
    ax.set_title("Lck Abundance", font = "Arial", fontsize = 16, fontweight = "bold")
    ax.set_ylabel("Corrected Lck Signal (a.u.)", font = "Arial", fontsize = 14, fontweight = "bold")
    ax.set_yticks([i*1000000 for i in range(8)])
    ax.set_yticklabels([f"{i}M" for i in range(8)], font = "Arial", fontsize = 12, fontweight = "bold")
    ax.set_xticks(centers)
    ax.set_xticklabels([samples[0][0], samples[1][0], samples[2][0]], font = "Arial", fontsize = 12, fontweight = "bold")#, rotation = 0, ha = "center", rotation_mode = "anchor")
    plt.savefig("lck_abundance.pdf", bbox_inches = "tight")
    plt.show()
    return None

#
#
#################################################################################################

    
if __name__ == "__main__":
    main()