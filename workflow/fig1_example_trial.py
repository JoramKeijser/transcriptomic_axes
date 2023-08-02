"""
Plot activity
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
import seaborn as sns
from sklearn.decomposition import PCA
from scipy.ndimage import convolve1d, gaussian_filter1d
sns.set_palette("colorblind")
sns.set_context("poster")


# Figure size
FONTSIZE = 7
FONTSIZE_SMALL = 7
LINEWIDTH = 2
width, height = 8, 8
INCH_TO_CM = 2.54


# Load activity from this session

# Load activity from this session
exp_folder = snakemake.input.experiment_dir + "/"
recording = exp_folder + "/Blank"
uids = np.squeeze(np.load(exp_folder + "neuron.UniqueID.npy"))
activity = np.load(recording + "/01/" + "frame.neuralActivity.npy")
states = np.squeeze(np.load(recording + "/01/" + "frame.states.npy"))
states_times = np.squeeze(np.load(recording + "/01/" + "frame.times.npy"))

prob_ttypes = np.squeeze(np.load(exp_folder + "neuron.prob_ttypes.npy"))
ttypes = np.loadtxt(exp_folder + "neuron.ttype.txt", dtype=str)

# Only keep cells Bugeon et al. could match to GABAergic Allen subtype,
# and with most likely cell type > 0.5
matched_to_allen = ~np.isnan(uids)
certain_assignment = np.max(prob_ttypes, 1) > 0.5
keep = matched_to_allen & certain_assignment
ecs = ttypes == 'EC'
ec_activity = activity[:, ecs]
gaba_activity = activity[:, keep]
# Order for plotting
running = gaba_activity[states==0].mean(0)
sync = gaba_activity[states==2].mean(0)
state_modulation = (running - sync) / (running + sync)
gaba_ttypes = ttypes[keep]
np.unique(gaba_ttypes)
subclasses = np.array([ttype.split("-")[0] for ttype in gaba_ttypes])
subclasses[subclasses=='Serpinf1'] = 'Vip'
subclass_order = ['Pvalb', 'Sst', 'Lamp5',  'Vip', 'Sncg']
plot_order = []
plot_order_subclass = {}
plot_order_subclass_alpha = {}
for subclass in subclass_order:
    subclass_idx = np.where(subclasses == subclass)[0]
    order = subclass_idx[np.argsort(state_modulation[subclass_idx])[::-1]]
    plot_order += list(order)
    plot_order_subclass[subclass] = order
    plot_order_subclass_alpha[subclass] = subclass_idx[np.argsort(gaba_ttypes[subclass_idx])[::-1]]
plt.plot(state_modulation[plot_order])

# Running times
running_speed = np.squeeze(np.load(recording + "/01/" + "running.speed.npy"))
running_times = np.squeeze(np.load(recording + "/01/" + "running.times.npy"))
eye_times = np.squeeze(np.load(recording + "/01/" + "eye.times.npy"))

activity_kernel_size = int(1 / (running_times[1] - running_times[0]))
activity_kernel = np.ones(activity_kernel_size) / activity_kernel_size
X_smooth = convolve1d(ec_activity, activity_kernel, axis=1)#, mode='same')
running_kernel_size = int(2 / (running_times[1] - running_times[0]))
running_kernel = np.ones(running_kernel_size) / running_kernel_size

# Where to start and end?
to_plot = states_times[states==2][np.where(np.diff(np.where(states == 2)[0])  > 1)[0]]
end_of_first = to_plot[0]
start_of_first = end_of_first - 40
bookends = {0: [100, 140], 1: [620, 660], 2: [start_of_first, end_of_first]}

# plotting params
dt_ratio = (eye_times[1] - eye_times[0] ) / (running_times[1]- running_times[0])
SIGMA = 5 # smoothing

# PCA on excitatory cells
X = ec_activity.T # cells x time
X_z = (X - X.mean(1, keepdims=True)) / X.std(1, keepdims=True)
X_smooth = gaussian_filter1d(X_z, sigma=2, axis=1)

pca = PCA(n_components=50)
pca.fit(ec_activity) # time x cells
ec_sort = np.argsort(pca.components_[0])

# Plot maps
colors = sns.color_palette("colorblind", 5)
cmaps = [LinearSegmentedColormap.from_list("",
                                            ["white", colors[i]]) for i in range(5)]
colors = sns.color_palette()
phase_names = {0:'Running', 2: 'Synchronized'}
bookends2 = {0: [100, 140], 2: [952.643, 992.643]}
n_cells = [len(plot_order_subclass[subclass]) for subclass in subclass_order]

n_cells = [n_cells[1], n_cells[1]] + n_cells
n_cells = np.array(n_cells)
height_ratios = n_cells / sum(n_cells)


gs = GridSpec(7, 2, height_ratios=height_ratios, hspace=0.1, wspace=0.1, left=0.25, right=0.9)
fig = plt.figure(figsize=(width/INCH_TO_CM, height/INCH_TO_CM), dpi=300)

for c, phase in enumerate([0, 2]):
    if c == 1:
        ax0 = fig.add_subplot(gs[0, c], sharey=ax0)
        ax1 = fig.add_subplot(gs[1, c], sharey=ax1)
    else:
        ax0 = fig.add_subplot(gs[0, c])
        ax1 = fig.add_subplot(gs[1, c])
    ax0.set_title(phase_names[phase], fontsize=FONTSIZE)

    # Running speed
    i0 = np.where(running_times >= bookends[phase][0]-5)[0][0]
    i1 = np.where(running_times >= bookends[phase][1]-5)[0][0]
    if c == 0:
        ax0.text(running_times[i0]-30,-.3, "Speed",
        fontsize=FONTSIZE, color='black', alpha=0.7)
        ax1.text(running_times[i0]-30,-.3, "Oscillation",
        fontsize=FONTSIZE, color='black', alpha=0.7)
        ax0.text(running_times[i0]-30,30, "State",
        fontsize=FONTSIZE, color='black')
    ax0.plot(running_times[i0:i1],
            gaussian_filter1d(running_speed[i0:i1],SIGMA*dt_ratio),
            lw=1, color='black', alpha=0.6)
    # Activity
    i0 = np.where(states_times >= bookends[phase][0]-5)[0][0]
    i1 = np.where(states_times >= bookends[phase][1]-5)[0][0]
    ax1.plot(states_times[i0:i1],
            X_smooth.T[i0:i1, ec_sort[:40]].mean(1),
            lw=1, color='black', alpha=0.6)

    for ax in [ax0, ax1]:
        ax.set_xticks([])
        ax.set_yticks([])


for r, subclass in enumerate(subclass_order):
    idx = plot_order_subclass_alpha[subclass]
    for c, phase in enumerate([2,0]):
        ax = fig.add_subplot(gs[r+2, c])
        if c==0:
            ax.set_ylabel(subclass,rotation=0,
                          fontsize=FONTSIZE, labelpad=20, color=colors[r])

        i0 = np.where(states_times >= bookends2[phase][0]-5)[0][0]
        i1 = np.where(states_times >= bookends2[phase][1]-5)[0][0]
        idx = plot_order_subclass_alpha[subclass]
        ax.pcolor(gaba_activity[i0:i1,idx].T, vmin = 0, vmax=450, cmap=cmaps[r])

        ax.set_yticks([])
        ax.set_xticks([])
        if (r==4) & (c == 0):
            Y0 = 0
            ax.hlines(-10, Y0, Y0+10 / np.diff(states_times)[0], lw=LINEWIDTH, color='black')
            ax.text(1, -35,  "10 s", fontsize=FONTSIZE_SMALL)
        if (r==0) & (c==1):
            x = i1-i0-1
            DY = 25
            ax.vlines(x, len(idx)-DY, len(idx)-DY+10, lw=LINEWIDTH, color='black')
            ax.text(x+15, len(idx)-37, "10 neurons",
                    fontsize=FONTSIZE_SMALL, rotation=90,  horizontalalignment='center')

sns.despine(bottom=True, left=True)
plt.savefig(snakemake.output.figure, dpi=300)
