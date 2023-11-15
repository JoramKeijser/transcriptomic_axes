"""
Simulate species-specific networks
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import array_to_latex as a2l
from src import constants, network_model
from src.network_model import simulate_variable_input

sns.set_context("poster")
cblind = sns.color_palette("colorblind")
palette = ["#666666", "#CC4125", cblind[0], cblind[1], cblind[2]]
sns.set_palette(palette=palette)


T0 = 300  # dont show init transient
# Shared across simulations
params = network_model.make_weights()
# Print weights as latex matrix

a2l.to_ltx(params["weights"], frmt="{:6.3f}", arraytype="pmatrix")
TIMESTEPS = 11000
time = np.arange(TIMESTEPS)
# Baseline rates. Sst active but not T0o active (Dend go)
baseline_rates = np.array([1, 1, 8, 4, 3])
x_base = network_model.find_inputs(params["weights"], baseline_rates)

# External and modulaT0ry inputs
x_var = np.zeros((TIMESTEPS, 5))
x_var[:, 0] = np.sin(time / 70) * 0.5 + 1
x_var[:, 1] = np.sin(time / 300) * 0.5 + 1

rates_dict = {}
# Mice
# ACh
Tm0, Tmend = 3300, 6700  # mod time start and end
x_mod = np.zeros((5,))
x_var[Tm0:Tmend, 2] = constants.PV_MOUSE_MOD
x_var[Tm0:Tmend, 3] = constants.SST_MOUSE_MOD
x_var[Tm0:Tmend, 4] = constants.VIP_MOUSE_MOD
rates_dict["mouse"] = simulate_variable_input(params, x_base, x_mod, x_var, dt=0.1)
# simulate_and_visualize(params, x_var, x_base, "mouse")

# Humans: only difference is no Pv inhibition
x_var[Tm0:Tmend, 2] = constants.PV_HUMAN_MOD
x_var[Tm0:Tmend, 3] = constants.SST_HUMAN_MOD
x_var[Tm0:Tmend, 4] = constants.VIP_HUMAN_MOD
rates_dict["human"] = simulate_variable_input(params, x_base, x_mod, x_var, dt=0.1)

# Turtle
x_var[Tm0:Tmend, 2] = constants.PV_TURTLE_MOD
x_var[Tm0:Tmend, 3] = constants.SST_TURTLE_MOD
x_var[Tm0:Tmend, 4] = constants.VIP_TURTLE_MOD
x_var[:Tm0, 3] = 5  # Sst baseline
x_var[Tmend:, 3] = 5
rates_dict["turtle"] = simulate_variable_input(params, x_base, x_mod, x_var, dt=0.1)

# Plot em
FONTSIZE = 15
ALPHA = 0.04
color = 'gray'
for species, figname in snakemake.output.items():
    rates = rates_dict[species]
    sns.set_palette("colorblind")
    fig, ax = plt.subplots(1, 3, figsize=(12, 3), sharex=True)
    for i in range(2):
        ax[0].plot(
            2 * (x_var[T0:, i] - x_var[T0:, i].mean(0)), alpha=0.7, lw=2, c=palette[i]
        )
        if species == "turtle":
            ax[1].fill_between(np.arange(0, Tm0)- T0, 0, 4, alpha=ALPHA, color=color)
            ax[2].fill_between(np.arange(0, Tm0)- T0, 1, 12, alpha=ALPHA, color=color)
            ax[1].fill_between(np.arange(Tmend,TIMESTEPS)- T0, 0, 4, alpha=ALPHA, color=color)
            ax[2].fill_between(np.arange(Tmend,TIMESTEPS)- T0, 1, 12, alpha=ALPHA, color=color)
        else:    
            ax[1].fill_between(np.arange(Tm0, Tmend)- T0, 0, 3.9, alpha=ALPHA, color=color)
            ax[2].fill_between(np.arange(Tm0, Tmend)- T0, 1, 12, alpha=ALPHA, color=color)
            #Tm0:Tmend

    ax[0].set_ylabel("Inputs (a.u.)", fontsize=FONTSIZE)
    ax[1].plot(rates[T0:, 0], lw=2, c=palette[0])
    ax[2].plot(rates[T0:, 2:], lw=2)
    ax[1].set_ylim([1, 4])
    ax[2].set_ylim([2, 12])
    ax[1].set_ylabel("PC rate (1/s)", fontsize=FONTSIZE)
    ax[2].set_ylabel("IN rate (1/s)", fontsize=FONTSIZE)
    for i in range(3):
        ax[i].set_xticks([0, 10000])
        ax[i].set_xticklabels([0, 10])
        ax[i].set_xlabel("Time (s)", fontsize=FONTSIZE)
    for i, (y, lw) in enumerate(zip([1.1, 3.9, 12], [3, 3, 6])):
        ax[i].hlines(y, Tm0 - T0, Tmend - T0, color="#81D2C7", lw=lw)
        # if i == 0:
        ax[i].text(Tm0 * 1.1, y * 1.1, "ACh", color="#81D2C7", fontsize=FONTSIZE)
        ax[i].tick_params(axis="both", which="major", labelsize=FONTSIZE)
    sns.despine()
    fig.tight_layout()
    plt.savefig(figname, dpi=300)
