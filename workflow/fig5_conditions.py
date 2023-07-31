"""
Visualize conditions for cholinergic modulation
like that in the mouse
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from src import network_model, constants
sns.set_context("poster")
sns.set_palette("colorblind")

params = network_model.make_weights()
VIP_MOD = 8 # fix Vip modulation
weights = params['weights']
n = weights.shape[0]
A = np.linalg.inv((np.eye(n) - weights))
xs = np.array([-0.5, 0.5])
mp_vals = np.array([-8, 8])

# Human: when is Delta r_pv < 0?
App, Aps, Apv = A[2, 2:] # Effective weights onto Pvalb neurons
# Line and arrow
plt.plot(mp_vals, -1/Aps * (App * mp_vals + Apv * VIP_MOD), lw=3)
plt.fill_between(mp_vals, 8, -1/Aps * (App * mp_vals + Apv * VIP_MOD), alpha=0.2)
plt.text(-6, 4, r"Pvalb inhibited", fontsize=20, color = sns.color_palette("colorblind")[0])

# Bla
Asp, Ass, Asv = A[3, 2:]
plt.plot(mp_vals, -1/Ass * (Asp * mp_vals + Asv * VIP_MOD), color=sns.color_palette()[1],
         label = 'dSst = 0', lw=3)
plt.fill_between(mp_vals, -8,
                 -1/Ass * (Asp * mp_vals + Asv * VIP_MOD), color=sns.color_palette()[1],
                 label = 'dSst = 0', lw=3, alpha = .2)
plt.xlabel("Pvalb  AChr expression")
plt.ylabel("Sst AChr expression")
plt.text(-1, -5, r"Sst inhibited", fontsize=20, color=sns.color_palette()[1])

# Annotate
plt.xlabel("Pvalb  AChr expression")
plt.ylabel("Sst AChr expression")
plt.ylim([-8, 8])
plt.xlim([-8, 8])
sns.despine()
plt.xticks([-8, 0, 8])
plt.yticks([-8, 0, 8])
plt.scatter(constants.PV_MOUSE_MOD, constants.SST_MOUSE_MOD,
            marker='*', color='black')
plt.tight_layout()
plt.savefig(snakemake.output.mouse, dpi=300)

# Human
plt.figure()
plt.plot(mp_vals, -1/Aps * (App * mp_vals + Apv * VIP_MOD), lw=3)
plt.plot(mp_vals, -1/Ass * (Asp * mp_vals + Asv * VIP_MOD),
         color=sns.color_palette()[1], label = 'dSst = 0', lw=3)
plt.xlabel("Pvalb  AChr expression")
plt.ylabel("Sst AChr expression")

# Annotate
plt.xlabel("Pvalb  AChr expression")
plt.ylabel("Sst AChr expression")
plt.ylim([-5, 5])
sns.despine()
plt.xticks([-8, 0, 8])
plt.yticks([-8, 0, 8])
plt.fill_between([-2, 2], 0, 8, color='gray', alpha=0.2)
plt.scatter(constants.PV_HUMAN_MOD, constants.SST_HUMAN_MOD, marker='*', color='black')
plt.tight_layout()
plt.savefig(snakemake.output.human, dpi=300)

# Turtle: when is Delta r_s > 0 ?
plt.plot(mp_vals, -1/Aps * (App * mp_vals + Apv * VIP_MOD), lw=3)
plt.plot(mp_vals, -1/Ass * (Asp * mp_vals + Asv * VIP_MOD),
         color=sns.color_palette()[1], label = 'dSst = 0', lw=3)
plt.xlabel("Pvalb  AChr expression")
plt.ylabel("Sst AChr expression")

# Annotate
plt.figure()
plt.plot(mp_vals, -1/Aps * (App * mp_vals + Apv * VIP_MOD), lw=3)
plt.plot(mp_vals, -1/Ass * (Asp * mp_vals + Asv * VIP_MOD),
          color=sns.color_palette()[1], label = 'dSst = 0', lw=3)
plt.xlabel("Pvalb  AChr expression")
plt.ylabel("Sst AChr expression")
plt.ylim([-5, 5])
sns.despine()
plt.xticks([-8, 0, 8])
plt.yticks([-8, 0, 8])
plt.fill_between([-2, 2], 0, -4, color='gray', alpha=0.2)
plt.scatter(constants.PV_TURTLE_MOD, constants.SST_TURTLE_MOD, marker='*', color='black')
plt.tight_layout()
plt.savefig(snakemake.output.turtle, dpi=300)
