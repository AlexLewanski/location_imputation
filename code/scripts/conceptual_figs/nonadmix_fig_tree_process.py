########################
### LOADING PACKAGES ###
########################

import tskit
import numpy as np
import random
import pandas as pd
 


###############################################
### LOADING AND INITIAL PROCESSING OF TREES ###
###############################################

### PROCESSING VARIABLES ###

ts_upload_path = '/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output/conceptual_fig/'
output_file_path = '/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output/conceptual_fig/'
output_tree_path = '/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output/conceptual_fig/'



### STEP 1: INITIAL STEPS
# 1A: LOAD INITIAL TREE SEQUENCE
# 1B: RECAPITATE TREE SEQUENCE
# ***NO SIMPLIFICATION

#1A
spatial_ts1 = tskit.load(ts_upload_path + "plot_sim_nonadmix.trees")


coalesce_info = [TREE.has_single_root for TREE in spatial_ts1.trees()]
     
spatial_ts1.num_trees
sum([not TREE.has_single_root for TREE in spatial_ts1.trees()])

#spatial_ts1.num_trees

#1B
#spatial_ts1 = pyslim.recapitate(
#    spatial_ts,
#    recombination_rate = 1e-8,
#    ancestral_Ne=500,
#    random_seed=921232)


random.seed(93317)

#sample indiviudals
sample_indiv_array = np.unique(spatial_ts1.nodes_individual[np.where(spatial_ts1.nodes_flags == 1)]) #sample individuals

# 2B: individuals in original population
subsample_sample = np.random.choice(sample_indiv_array, size=12, replace=False)

subsample_sample = np.random.choice(sample_indiv_array, size=12, replace=False)

subsamp_node_list = [spatial_ts1.individual(x).nodes[0] for x in subsample_sample]



#2C: simplify
#spatial_ts1_reduced = spatial_ts1.simplify(np.concatenate(subsamp_node_list),
#                                           keep_unary = False)
spatial_ts1_reduced = spatial_ts1.simplify(subsamp_node_list,
                                           keep_unary = False)

#spatial_ts1_reduced.num_trees
#spatial_ts1_reduced.keep_intervals(intervals = np.array([366965.0, 366966.0]),
#                                   record_provenance=True,
#                                   simplify = False)


#tree_inds = [1, 13, 21, 46, 63]
tree_inds = [1, 6, 28, 51, 66]

node_rows = []
edge_rows = []

for IND in tree_inds:
    tree1 = spatial_ts1_reduced.at_index(IND)
    
    for u in tree1.nodes():
        #print(u)
        n = spatial_ts1_reduced.node(u)
        print(n)
        node_rows.append({
            "tree_ind": IND,
            "node_id": u,
            "time": n.time,
            "x": spatial_ts1_reduced.individual(spatial_ts1_reduced.node(u).individual).location[0],
            "y": spatial_ts1_reduced.individual(spatial_ts1_reduced.node(u).individual).location[1]
        })

    for u in tree1.nodes():
        parent = tree1.parent(u)
        if parent != tskit.NULL:
            edge_rows.append({
                "tree_ind": IND,
                "parent": parent,
                "child": u
            })

nodes_df = pd.DataFrame(node_rows)
edges_df = pd.DataFrame(edge_rows)

span_info = []
for i,TREE in enumerate(spatial_ts1_reduced.trees()):
    if i in tree_inds:
        tree_group = 't' + str(i)
    else:
        tree_group = 'n'
    span_info.append({
        "left": TREE.interval.left,
        "right": TREE.interval.right,
        "group": tree_group
    })

span_df = pd.DataFrame(span_info)


nodes_df.to_csv("/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output/conceptual_fig/nodes_nonadmix.csv", index=False)
edges_df.to_csv("/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output/conceptual_fig/edges_nonadmix.csv", index=False)
span_df.to_csv("/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output/conceptual_fig/span_info_nonadmix.csv", index=False)
