import numpy as np
import anndata as ad

def load_activity(base_dir, exp_dir, recording_dir):
    """
    Load activity data from Bugeon et al. 
    Arguments:
        base_dir: directory with all the data
        exp_dir: directory of particular experiment 
        recording_dir: one recording

    Returns:
        annotated data frame (neurons x time)
    """
    # Load the neuron data
    uids = np.squeeze(np.load(exp_dir + "neuron.UniqueID.npy"))
    counts = np.load(exp_dir + "neuron.gene_count.npy")
    prob_ttypes = np.squeeze(np.load(exp_dir + "neuron.prob_ttypes.npy"))
    ttypes = np.loadtxt(exp_dir + "neuron.ttype.txt", dtype=str)

    # Load session data
    activity = np.load(recording_dir + "/01/" + "frame.neuralActivity.npy")
    states = np.squeeze(np.load(recording_dir + "/01/" + "frame.states.npy"))
    states_times = np.squeeze(np.load(recording_dir + "/01/" + "frame.times.npy"))
    running_speed = np.squeeze(np.load(recording_dir + "/01/" + "running.speed.npy"))
    running_times = np.squeeze(np.load(recording_dir + "/01/" + "running.times.npy"))

    # Subset interneurons
    matched_to_allen = ~np.isnan(uids)
    certain_assignment = np.max(prob_ttypes, 1) > 0.5
    keep = matched_to_allen & certain_assignment 
    gaba_activity = activity[:, keep]
    gaba_ttypes = ttypes[keep]
    states_names = np.loadtxt(base_dir + "states.names.txt", dtype=str, delimiter=',')
    
    # Put in anndata
    adata = ad.AnnData(X=gaba_activity)
    adata.obs['states'] = states_names[np.array(states, dtype=int)]
    adata.obs['states'] = adata.obs['states'].astype("category")
    adata.var['Subtype'] = gaba_ttypes
    adata.var['Subclass'] = [ttype.split("-")[0] for ttype in gaba_ttypes]
    adata.var['Subclass'] = adata.var['Subclass'].astype("category")
    # Put Serpinf1 neurons into Vip subclass
    if 'Serpinf1' in adata.var['Subclass'].cat.categories:
        adata.var.loc[adata.var['Subclass'] == 'Serpinf1', 'Subclass'] = 'Vip'
        adata.var['Subclass'] = adata.var['Subclass'].cat.remove_categories(['Serpinf1'])
    # Add missing categories
    subclasses = ['Pvalb', 'Sst', 'Lamp5', 'Vip', 'Sncg']
    missing_categories = [subclass for subclass in subclasses if subclass not in np.unique(adata.var['Subclass'])]
    adata.var['Subclass'] = adata.var['Subclass'].cat.add_categories(missing_categories)
    adata.var['Subclass'] = adata.var['Subclass'].cat.reorder_categories(subclasses)
    adata.obs['time'] = states_times
    adata.uns['Running speed'] = running_speed
    adata.uns['time'] = running_times

    # Neurons = variables
    adata.var.index = np.array(uids[keep], dtype=int)
    adata.var.index.rename("uid", inplace=True)

    # Average activity in state
    for state in states_names:
        adata.var[state] = adata.X[adata.obs['states'] == state].mean(0)

    # Difference in avg activity in different states
    adata.var['State modulation'] = (adata.var['Running'] - adata.var['Stationary Synchronized']) /  \
                                (adata.var['Running'] + adata.var['Stationary Synchronized'])
    adata.var['Running modulation'] = (adata.var['Running'] - adata.var['Stationary Desynchronized']) /  \
                                (adata.var['Running'] + adata.var['Stationary Desynchronized'])
    adata.var['Stationary modulation'] = (adata.var['Stationary Desynchronized'] - adata.var['Stationary Synchronized']) /  \
                                (adata.var['Stationary Desynchronized'] + adata.var['Stationary Synchronized'])

    return adata

def organize_subclass_labels(adata):
    """
    Add labels for missing subclasses, 
    order according to mouse tPC1
    """
    subclass_order = ['Pvalb', 'Sst', 'Lamp5', 'Vip', 'Sncg', 'Meis2']
    adata.obs['Subclass'] = adata.obs['Subclass'].astype("category")
    missing_subclasses = [subclass for subclass in subclass_order if subclass not in adata.obs['Subclass'].cat.categories]
    adata.obs['Subclass'] = adata.obs['Subclass'].cat.add_categories(missing_subclasses)
    adata.obs['Subclass'] = adata.obs['Subclass'].cat.reorder_categories(subclass_order)
    return adata