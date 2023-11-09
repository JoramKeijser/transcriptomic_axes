"""
Schematics of cross-variance 
"""

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from sklearn.decomposition import PCA
import os
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_palette("colorblind")
sns.set_context("poster")
sc.settings.dpi_save= 300

# Mouse vs human PC
n = 2000
fontsize = 30
cov = np.array([[1, 0.6], [0.6, .5]]) * 0.2
evecs = np.linalg.eig(cov)[1]
evecs *= 1.5 #
rng = np.random.RandomState(123)
sns.set_palette("Set2")
x = rng.multivariate_normal(np.zeros((2, )), cov, size = (n, ))
# Mouse
plt.scatter(x[:,0], x[:,1], alpha=0.2, s=20, color='gray')
plt.plot([0, evecs[0,0]], [0,evecs[1,0]], color=sns.color_palette()[0], lw=6)
plt.text(evecs[0,0]+0.1, evecs[1,0]+0.1, "Mouse \n tPC1", fontsize=fontsize, color=sns.color_palette()[0],fontweight='bold')

# Rotate human PCs
theta = 20 * np.pi / 180
R = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
evecs2 = R@evecs
evecs2  *= 0.4
plt.plot([0, evecs2[0,0]], [0,evecs2[1,0]], color=sns.color_palette()[1], lw=6)
plt.text(evecs2[0,0]-0.7, evecs2[1,0]+0.2, "Human \n tPC1", fontsize=fontsize, color=sns.color_palette()[1],fontweight='bold')


sns.despine(bottom=True, left=True)
plt.xticks([])
plt.yticks([])
plt.tight_layout()
plt.savefig(snakemake.output.pc_sketch) #"../figures/figure2/pcs_sketch"


# 2. Predict mouse modulation from human PC
# State modulation: gradient
gradient = x@evecs[:,0] 
noise = rng.normal(size=n) * 0.4
noisy_gradient = gradient + noise
plt.scatter(x[:,0], x[:,1], alpha=0.2, s=20,c=noisy_gradient)
plt.scatter(np.linspace(0, evecs[0,0],n), np.linspace(0,evecs[1,0], n), c=np.linspace(0, 1, n), 
            s=10, cmap="viridis")
plt.text(evecs[0,0]+0.1, evecs[1,0]-0.2, "State \n mod", fontsize=fontsize, color=sns.color_palette()[0],fontweight='bold')

# Human PC, as before
plt.plot([0, evecs2[0,0]], [0,evecs2[1,0]], color=sns.color_palette()[1], lw=6)
plt.text(evecs2[0,0]-1, evecs2[1,0]+0.2, "Human \n tPC1", fontsize=fontsize, color=sns.color_palette()[1],fontweight='bold')

sns.despine(bottom=True, left=True)
plt.xticks([])
plt.yticks([])
plt.tight_layout()
plt.savefig(snakemake.output.pc_gradient) #"../figures/figure2/pcs_gradient"