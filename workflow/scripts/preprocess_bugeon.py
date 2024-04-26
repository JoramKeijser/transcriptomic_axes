"""
Extract the transcriptomic and activity data
from the raw Bugeon et al. files downloaded from
https://doi.org/10.6084/m9.figshare.19448531.v1
Save as annotated data frame to /results/
"""
import os
import numpy as np
import seaborn as sns
import pandas as pd
import anndata as ad

sns.set_palette("colorblind")
sns.set_context("poster")


base_dir = snakemake.input.data_dir + "/"
subdirs = [d for d in os.listdir(base_dir) if os.path.isdir(base_dir + d)]

# Dictionaries to hold cell-specific data, indexed by cell (unique identifier)
# activity data for each behavioral state
# coded as (0,1,2) = (Running, Stat Desynch., Stat Synch)
all_running = {}
all_desync = {}
all_sync = {}
# time spend in each behavioral state
all_running_pct = {}
all_sync_pct = {}
all_desync_pct = {}
# cell type
all_ttypes = {}
all_maxprobs = {}
# RNA counts
all_counts = {}

# Loop over sub folders
for subdir in subdirs:
    full_subdir = base_dir + subdir + "/"
    experiments = os.listdir(full_subdir)
    for experiment in experiments:
        exp_folder = full_subdir + experiment + "/"

        # Load activity data
        recording_folders = [f.path for f in os.scandir(exp_folder) if f.is_dir()]
        for recording in recording_folders:
            uids = np.squeeze(np.load(exp_folder + "neuron.UniqueID.npy"))
            counts = np.load(exp_folder + "neuron.gene_count.npy")
            prob_ttypes = np.squeeze(np.load(exp_folder + "neuron.prob_ttypes.npy"))
            ttypes = np.loadtxt(exp_folder + "neuron.ttype.txt", dtype=str)
            subfolder = (
                recording + "/01/"
            )  # some also have /02/ but the cell info doesn't match
            activity = np.load(subfolder + "/frame.neuralActivity.npy")
            states = np.squeeze(np.load(subfolder + "/frame.states.npy"))
            # Only keep cells Bugeon et al. could match to GABAergic Allen subtype,
            # with max posterior cell type > 0.5
            matched_to_allen = ~np.isnan(uids)
            certain_assignment = np.max(prob_ttypes, 1) > 0.5
            keep = matched_to_allen & certain_assignment
            # subset cells, and select imaging data from each behavioral state
            counts = counts[keep]
            running = activity[states == 0][:, keep]
            desync = activity[states == 1][:, keep]
            sync = activity[states == 2][:, keep]
            ttypes = ttypes[keep]
            prob_ttypes = prob_ttypes[keep]
            uids = uids[keep]
            # Compute percentage in each state - see states.names.txt
            running_pct = np.mean(states == 0) * 100
            desync_pct = np.mean(states == 1) * 100
            sync_pct = np.mean(states == 2) * 100

            # Loop over cells (unique identifiers, uids)
            for i, uid in enumerate(uids):
                # Can just over-write t-info, it's the same for all sessions
                all_maxprobs[uid] = np.max(prob_ttypes[i])
                all_ttypes[uid] = ttypes[i]
                all_counts[uid] = counts[i]
                # Now check if we should use activity data for this session
                # because we didn't see the neuron before, or if this session has more running
                # than the existing one
                neuron_exists = uid in all_running.keys()
                if not (neuron_exists) or sync_pct >= all_sync_pct[uid]:
                    # % state per neuron
                    all_running_pct[uid] = running_pct
                    all_desync_pct[uid] = desync_pct
                    all_sync_pct[uid] = sync_pct
                    # Activity
                    all_running[uid] = running[:, i]
                    all_desync[uid] = desync[:, i]
                    if sync.shape[0] > 0:
                        all_sync[uid] = sync[:, i]
                    else:
                        all_sync[uid] = np.empty((1,))
                        all_sync[uid].fill(np.nan)

# Dict to array
running_pct = np.array(list(all_running_pct.values()))
sync_pct = np.array(list(all_sync_pct.values()))
desync_pct = np.array(list(all_desync_pct.values()))
print(
    "% of time in running, stat. sync, and stat. desync state: \
      {running_pct.mean():1.0f}, {sync_pct.mean():1.0f}, {desync_pct.mean():1.0f}"
)
n_cells = len(running_pct)
print(f"Found {n_cells} GABAergic cells")
# Make sure percentages add up to 100
assert np.allclose(sync_pct + running_pct + desync_pct, 100)

# Compute average activity during each state
uids = np.array(list(all_running_pct.keys()))
ttypes = np.array(list(all_ttypes.values()))
avg_running = np.array([activity.mean() for activity in all_running.values()])
avg_desync = np.array([activity.mean() for activity in all_desync.values()])
avg_sync = np.array([activity.mean() for activity in all_sync.values()])
# State modulation: normalized activity difference
mod = (avg_running - avg_desync) / (avg_running + avg_desync)
print(
    f"Found {np.isnan(avg_running).sum()} cells not recorded during running state"
)  # 0
print(
    f"Found {np.isnan(avg_sync).sum()} cells not recorded during synchronized state"
)  # 139
print(
    f"Found {np.isnan(avg_desync).sum()} cells not recorded during desynchronized state"
)  # 0

# Put all metadata into data frame
df = pd.DataFrame(np.array([uids, mod]).T, columns=["uid", "state_modulation"])
df["Subtype"] = ttypes
df["uid"] = df["uid"].astype(int)
df["Subclass"] = [ttype.split("-")[0] for ttype in ttypes]

# Final cell filter
sparse_types = df["Subtype"].value_counts()[df["Subtype"].value_counts() < 3].index
print(
    "Exclude",
    sum(df["Subtype"].isin(sparse_types)),
    "neurons of subtype with < 3 cells",
)  # 8, as listed in paper
df_exclude = df[df["Subtype"].isin(sparse_types)]
df = df[~df["Subtype"].isin(sparse_types)]
# As categories
df = df.set_index(df["uid"])
df = df.drop(["uid"], axis=1)
# Reorder columns
df = df[["Subclass", "Subtype", "state_modulation"]]
# Put Serpinf1 into Vip subclass
df["Subclass"] = df["Subclass"].replace({"Serpinf1": "Vip"})
df["Subclass"] = df["Subclass"].astype("category")
subclass_order = ["Pvalb", "Sst", "Lamp5", "Vip", "Sncg"]
df["Subclass"] = df["Subclass"].cat.reorder_categories(subclass_order)
df["Subclass"] = df["Subclass"].astype("category")


#  RNA counts (already normalized) as data frame
counts = np.array(list(all_counts.values()))
gene_names = np.loadtxt(base_dir + "genes.names.txt", dtype=str)
counts = pd.DataFrame(counts, columns=gene_names)
counts.set_index(uids, inplace=True)
# also filter cells of the counts
counts = counts[counts.index.isin(df.index)]
assert counts.shape[0] == df.shape[0]  # now the number of cells should match

num_ttypes = len(df["Subtype"].value_counts())
num_cells = df.shape[0]
print(f"Left with {num_cells} of {num_ttypes} subtypes")  # 1065 and 35, as in paper


# Create annotated data frame
bugeon = ad.AnnData(counts.to_numpy(), obs=df, dtype=np.float64)
bugeon.var_names = counts.columns

# Order subclasses by state modulation
subclass_order = ["Pvalb", "Sst", "Lamp5", "Vip", "Sncg"]
bugeon.obs["Subclass"] = bugeon.obs["Subclass"].cat.reorder_categories(subclass_order)

# Npy is missing (NaN) for some cells
missing_cells = np.sum((np.isnan(bugeon.X) > 0).sum(1) > 0)
missing_genes = np.sum((np.isnan(bugeon.X) > 0).sum(0) > 0)
print(f"{missing_cells} cells with missing counts in {missing_genes} genes")
gene_idx = np.where((np.isnan(bugeon.X) > 0).sum(0) > 0)[0]
print(f"Genes with missing counts: {np.array(bugeon.var_names[gene_idx])}")
# Replace by subtype-specific median
for gene in bugeon.var_names[gene_idx]:
    for subtype in np.unique(bugeon.obs["Subtype"]):
        idx = (bugeon.obs["Subtype"] == subtype) & np.isnan(bugeon[:, gene].X[:, 0])
        if idx.sum() > 0:
            bugeon[idx, gene].X = np.nanmedian(
                bugeon[bugeon.obs["Subtype"] == subtype, gene].X
            )

print(f"Save anndata to {snakemake.output.anndata}")
# sort by uid
bugeon = bugeon[np.argsort(np.array(bugeon.obs_names, dtype=int))]
bugeon.write_h5ad(snakemake.output.anndata)
