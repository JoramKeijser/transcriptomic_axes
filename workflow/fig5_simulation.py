import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import array_to_latex as a2l
from src import constants
sns.set_context("poster")
cblind = sns.color_palette("colorblind")
palette = ['#666666', '#CC4125', cblind[0], cblind[1], cblind[2]]
sns.set_palette(palette=palette)
from src import network_model
FIGDIR = "./figures/figure6/"


def simulate_and_visualize(params, x_var, x_base, species, t0 = 300, save = False):
    x_mod = np.zeros((5, ))
    rates = network_model.simulate_variable_input(params, x_base, x_mod, x_var, dt = 0.1)
    
    if save:
        fig, ax = plt.subplots(2, 1, figsize=(6, 6), sharex=True)
        ax[0].plot(2*(x_var[t0:,:2] - x_var[t0:,:2].mean(0)), alpha = 0.7, lw=3)
        ax[1].plot(rates[t0:,0], lw=2)

        for i, y in enumerate([1.1, rates[t0:,0].max()]):
            ax[i].hlines(y, 3000, 8000,color = "#81D2C7")
            if i == 0:
                ax[i].text(4500, y*1.1, "ACh", color = "#81D2C7")

        ax[0].set_yticks([0])
        ax[0].set_ylabel("inputs (norm.)")
        ax[1].set_ylabel("PC rate")
        ax[1].set_xlabel("time (ms)")
        sns.despine()
        fig.tight_layout()
        plt.savefig(FIGDIR + f"simulation_{species}", dpi=300)
    return rates


def main():
    t0 = 300 # dont show init transient
    # Shared across simulations
    params = network_model.make_weights()
    # Print weights as latex matrix
    
    a2l.to_ltx(params['weights'], frmt = '{:6.3f}', arraytype = 'pmatrix')
    time_steps = 11000
    time = np.arange(time_steps)
    # Baseline rates. Sst active but not too active (Dend go)
    baseline_rates = np.array([1,1,8,4,3]) 
    x_base = network_model.find_inputs(params['weights'], baseline_rates)

    # External and modulatory inputs
    x_var = np.zeros((time_steps, 5))
    x_var[:,0] = np.sin(time / 70)*.5+1
    x_var[:,1] = np.sin(time / 300)*.5+1

    rates_dict = {}
    # Mice
    # ACh
    Tm0, Tmend = 3300, 6700 # mod time start and end
    x_var[Tm0:Tmend,2] = constants.PV_MOUSE_MOD 
    x_var[Tm0:Tmend,3] = constants.SST_MOUSE_MOD
    x_var[Tm0:Tmend,4] = constants.VIP_MOUSE_MOD
    rates_dict['mouse'] = simulate_and_visualize(params, x_var, x_base, "mouse")

    # Humans: only difference is no Pv inhibition
    x_var[Tm0:Tmend,2] = constants.PV_HUMAN_MOD 
    x_var[Tm0:Tmend,3] = constants.SST_HUMAN_MOD 
    x_var[Tm0:Tmend,4] = constants.VIP_HUMAN_MOD 
    rates_dict['human'] = simulate_and_visualize(params, x_var, x_base, "human")

    # Turtle
    x_var[Tm0:Tmend,2] = constants.PV_TURTLE_MOD 
    x_var[Tm0:Tmend,3] = constants.SST_TURTLE_MOD 
    x_var[Tm0:Tmend,4] = constants.VIP_TURTLE_MOD 
    x_var[:Tm0,3] = 5 # Sst baseline
    x_var[Tmend:,3] = 5
    rates_dict['turtle'] = simulate_and_visualize(params, x_var, x_base, "turtle")

    # Plot em
    fs = 15
    for species, figname in snakemake.output.items():
        print(species, figname)
        rates = rates_dict[species]
        sns.set_palette("colorblind")
        fig, ax = plt.subplots(1, 3, figsize=(12, 3), sharex=True)
        for i in range(2):
            ax[0].plot(2*(x_var[t0:,i] - x_var[t0:,i].mean(0)), alpha = 0.7, lw=2, c=palette[i])

        ax[0].set_ylabel("Inputs (a.u.)", fontsize=fs)
        ax[1].plot(rates[t0:,0], lw=2,  c = palette[0])
        ax[2].plot(rates[t0:,2:], lw=2)
        ax[1].set_ylim([1, 4])
        ax[2].set_ylim([2, 12])
        ax[1].set_ylabel("PC rate (1/s)", fontsize=fs)
        ax[2].set_ylabel("IN rate (1/s)", fontsize=fs)
        for i in range(3):
            ax[i].set_xticks([0, 10000])
            ax[i].set_xticklabels([0, 10])
            ax[i].set_xlabel("Time (s)", fontsize=fs)
        for i, (y, lw) in enumerate(zip([1.1, 3.9, 12], [3, 3, 6])):
            ax[i].hlines(y, Tm0-t0, Tmend-t0,color = "#81D2C7", lw = lw)
            #if i == 0:
            ax[i].text(Tm0*1.1, y*1.1, "ACh", color = "#81D2C7", fontsize=fs)
            ax[i].tick_params(axis='both', which='major', labelsize=fs)
        sns.despine()
        fig.tight_layout()
        plt.savefig(figname, dpi=300)

if __name__=="__main__":
    main()