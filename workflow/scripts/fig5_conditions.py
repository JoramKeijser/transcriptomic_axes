"""
Visualize conditions on channel densities
for cholinergic effect to be qualitatively 
consistent with the mouse data
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from src import network_model, constants
sns.set_context("poster")
sns.set_palette("colorblind")

def plot_conditions(species, VIP_MOD = 8, verbose = False):
    fig, ax = plt.subplots()
    sns.despine()
    mp_vals = np.arange(-4, 4, 0.1)
    pv_conditions = np.zeros((2, len(mp_vals)))
    sst_conditions = np.zeros((2, len(mp_vals)))

    for i, dendrite in enumerate(['on', 'off']):
        params = network_model.make_weights()
        weights = params["weights"]
        if dendrite == 'off':
            weights[0,1] = 0
        A = np.linalg.inv((np.eye(weights.shape[0]) - weights))
        App, Aps, Apv = A[2, 2:]  # Effective weights onto Pvalb neurons
        Asp, Ass, Asv = A[3, 2:]
        if verbose:
            print("Dendrite:", dendrite)
            print(f"A_pp: {A[2,2]:0.3f}, A_ps: {A[2,3]:0.3f}, A_pv: {A[2,4]:0.3f}")    
        pv_conditions[i] = -1 / Aps * (App * mp_vals + Apv * VIP_MOD)
        sst_conditions[i] = -1 / Ass * (Asp * mp_vals + Asv * VIP_MOD)

    # Pick most stringent condition at each point
    pv_conditions = np.max(pv_conditions, axis=0)
    sst_conditions = np.max(sst_conditions, axis=0)

    plt.plot(mp_vals, pv_conditions, color=sns.color_palette()[0], lw=2)
    plt.plot(
        mp_vals,
        sst_conditions,
        color=sns.color_palette()[1],
        label="dSst = 0",
        lw=2,
    )
    if species.lower() == "mouse":
        plt.scatter(constants.PV_MOUSE_MOD, constants.SST_MOUSE_MOD, marker="*", color="black")
        plt.text(-1, -5, r"Sst inhibited", fontsize=20, color=sns.color_palette()[1])
        plt.text(
            -1, 5.5, r"Pvalb inhibited", fontsize=20, color=sns.color_palette("colorblind")[0])

        # Shade for mouse
        plt.fill_between(mp_vals, 7, pv_conditions, 
                    alpha=0.2, 
                    color= sns.color_palette()[0])
        plt.fill_between(
            mp_vals,
            -8,
            sst_conditions,
            color=sns.color_palette()[1],
            label="dSst = 0",
            lw=3,
            alpha=0.2,
        )
    elif species.lower() == "human":
        plt.scatter(constants.PV_HUMAN_MOD, constants.SST_HUMAN_MOD, marker="*", color="black")
        plt.fill_between([-2, 2], 0, 6.8, color="gray", alpha=0.2, linestyle=":")

    elif species.lower() == "turtle":
        plt.fill_between([-2, 2], 0, -4, color="gray", alpha=0.2, linestyle = ":")
        plt.scatter(
            constants.PV_TURTLE_MOD, constants.SST_TURTLE_MOD, marker="*", color="black"
        )
    else:
        print(f"Species {species} not implemented")

    
    ax.set_xlim([-3, 3])
    ax.set_ylim([-7, 7])
    ax.set_xticks([-3, 0, 3])
    ax.set_yticks([-7, 0, 7])
    ax.set_xlabel("Pvalb  AChr expression")
    ax.set_ylabel("Sst AChr expression")
    fig.tight_layout()
    return fig, ax, mp_vals, pv_conditions, sst_conditions

# Fix vip modulation to be large, positive. 
# Different values give qualitatively similar results
VIP_MOD = 8 
fig, ax, mp_vals, pv_conditions, sst_conditions = plot_conditions("mouse", VIP_MOD, True)
fig.savefig(snakemake.output.mouse, dpi=300)
fig, ax, mp_vals, pv_conditions, sst_conditions = plot_conditions("human", VIP_MOD)
plt.savefig(snakemake.output.human, dpi=300)
fig, ax, mp_vals, pv_conditions, sst_conditions = plot_conditions("turtle", VIP_MOD)
plt.savefig(snakemake.output.turtle, dpi=300)