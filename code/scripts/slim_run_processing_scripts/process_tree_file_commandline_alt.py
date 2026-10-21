########################
### LOADING PACKAGES ###
########################

import argparse
import tskit
import msprime
import numpy as np
import pyslim
import random
import gzip
from pathlib import Path

#import sys
#import msprime
#import networkx as nx
#from IPython.display import SVG
#import math



########################
### CUSTOM FUNCTIONS ###
########################

#map from new to old node ids, taken from tskit website: https://tskit.dev/tutorials/advanced_simplification.html
def invert_map(node_mapping):
    kept = node_mapping != tskit.NULL
    indexes = node_mapping[kept]  # indexes are guaranteed 0..N-1
    rev_map = np.full_like(indexes, tskit.NULL)
    rev_map[indexes] = np.flatnonzero(kept)
    return rev_map

def extract_samples_multipop(ts, exclude_pops):
    sample_node_list = []
    sample_node_pop = []

    for node_id in ts.samples():
        population_id = ts.node(node_id).population
        if population_id not in exclude_pops:
            sample_node_list.append(node_id)
            sample_node_pop.append(population_id)

    return {
        "sample_node_list": sample_node_list,
        "pop_list": sample_node_pop,
    }

def extract_ped_founder_anc(ts, focal_node):
    ancestor_id_list = []
    for TREE in ts.trees():
        node = focal_node
        subpop = ts.node(node).population
        while subpop != 1:
            node_init = TREE.parent(node)
            subpop = ts.node(node_init).population

            #the nodes of the founder individuals moved to the new subpop will still have a subpop of 1, so we
            #are basically looking for the first node to have a subpop of one going back in time. This is slightly
            #different then if we looked up the subpop from the individual associated with the node, whose subpop
            #would apparently be the new subpop (I'm not 100% sure why). If we were looking up the subpop from the
            #individual rather than node, we would have to do this slightly differently.
            node = node_init
        ancestor_slim_id = ts.individual(ts.node(node).individual).metadata['pedigree_id']
        ancestor_id_list.append(ancestor_slim_id)
    return(ancestor_id_list)


def ancestry_intervals(ts, anc_list, focal_node = None):
    id_counter = 0
    pedig_id_list = []
    left_interv_list = []
    right_interv_list = []

    total_ind = len(anc_list) - 1
    for index,pedig_id in enumerate(anc_list):
        if id_counter == 0:
            pedig_id_list.append(pedig_id)
            left_interv_list.append(ts.at_index(index).interval.left)

        #the last tree
        if index == total_ind:
            right_interv_list.append(ts.at_index(index).interval.right)
            break

        if pedig_id != anc_list[index + 1]:
            right_interv_list.append(ts.at_index(index).interval.right)
            id_counter = 0
        else:
            id_counter += 1

    return({"focal_node": focal_node,
            "ancestor_list": pedig_id_list,
            "left_interv_list": left_interv_list,
            "right_interv_list": right_interv_list})



#############################
### PARSE ARGUMENT INPUTS ###
#############################

parser = argparse.ArgumentParser(description="Process tree sequence.")

# Positional arguments
parser.add_argument("--tree_sep", "-t", type=int, help="separation for tree sampling", required=True)
parser.add_argument("--sample_size", "-s", type=int, help="separation for tree sampling", required=True)
parser.add_argument("--ts_upload_path", "-u", help="upload path for processed tree sequence", required=True)
parser.add_argument("--ts_name", "-m", help="name of tree sequence to upload", required=True)
parser.add_argument("--output_file_path", "-o", help="output path for files", required=True)
parser.add_argument("--output_tree_path", "-p", help="output path for newick trees", required=True)
parser.add_argument(
    "--tree_output",
    choices=["separate", "gzipped", "both"],
    default="gzipped",
    help="Newick output format: separate files, one gzipped file, or both (default: gzipped)",
    required=True
)
parser.add_argument("--output_prefix", "-f", help="output name for newick trees", required=True)
parser.add_argument("--recomb_rate", "-r", type=float, help="recombination rate for msprime simulation", required=True)
parser.add_argument("--ne", "-n", type=float, help="effective size for msprime simulation", required=True)
parser.add_argument("--seed", "-d", type=int, help="random seed", required=True)



###################################
### PROCESS INPUT TREE SEQUENCE ###
###################################
def main(args):

    #initial checks
    if args.tree_sep <= 0:
        raise ValueError("--tree_sep must be a positive integer")

    #setting up seeds
    random.seed(args.seed)
    np.random.seed(args.seed)

    ###############################################
    ### LOADING AND INITIAL PROCESSING OF TREES ###
    ###############################################

    ### STEP 1: INITIAL STEPS
    # 1A: LOAD INITIAL TREE SEQUENCE
    # 1B: RECAPITATE TREE SEQUENCE
    # ***NO SIMPLIFICATION

    #1A
    #make path more robust to when or not input pah contains trailing back slash
    ts_path = Path(args.ts_upload_path) / args.ts_name
    spatial_ts = tskit.load(ts_path)
    #spatial_ts = tskit.load(args.ts_upload_path + args.ts_name)

    #1B
    spatial_ts1 = pyslim.recapitate(
        spatial_ts,
        recombination_rate = args.recomb_rate,
        ancestral_Ne = args.ne,
        random_seed=args.seed)

    ### STEP 2: INITIAL SIMPLIFICATION INCLUDING ADMIXED INDIVS + OTHERS
    # 2A: IDENTIFY ADMIXED INDIVS
    # 2B: IDENTIFY NON-ADMIXED INDIVS TO BE SAMPLED
    # 2C: SIMPLIFY BASED ON ADMIXED AND SAMPLED INDIVS (KEEP UNARY NODES!)

    #sample individuals
    sample_indiv_array = np.unique(spatial_ts1.nodes_individual[np.where(spatial_ts1.nodes_flags == 1)]) #sample individuals

    #2A
    #admix_nodes = spatial_ts1.samples(population=2) #nodes in admixed population
    admix_nodes_dict = extract_samples_multipop(ts = spatial_ts1, exclude_pops = [1])
    admix_nodes_original=admix_nodes_dict["sample_node_list"]
    admix_indivs = list(sorted(set([spatial_ts1.node(i).individual for i in admix_nodes_original]))) #individuals in the admixed pop


    #info on ancestry node
    anc_node_info_list = []
    for ADMIX_NODE in admix_nodes_original:
        anc_list = extract_ped_founder_anc(ts = spatial_ts1, focal_node = ADMIX_NODE)
        anc_node_info_list.append(ancestry_intervals(ts = spatial_ts1, anc_list = anc_list, focal_node = ADMIX_NODE))

    #individuals in original population
    non_admixed_indiv_array = sample_indiv_array[~np.isin(sample_indiv_array, admix_indivs)]

    # 2B: individuals in original population
    if len(non_admixed_indiv_array) < args.sample_size:
        raise ValueError("The number of available individuals is less than the requested sample size!")

    subsample_sample = np.random.choice(non_admixed_indiv_array, size=args.sample_size, replace=False)
    subsamp_node_list = [spatial_ts1.individual(x).nodes for x in subsample_sample]


    ### SIMPLIFY DOWN TO EACH FOCAL INDIVIDUAL
    simplify_tree_list = []

    for ADMIX_INDIV in admix_indivs:
        focal_nodes = spatial_ts1.individual(ADMIX_INDIV).nodes

        spatial_ts1_reduced,node_map = spatial_ts1.simplify(np.concatenate([np.concatenate(subsamp_node_list), focal_nodes]),
                                                            map_nodes = True,
                                                            filter_individuals = False)
        reverse_map = invert_map(node_map)

        #new admix node ids
        focal_nodes_postsimpl = extract_samples_multipop(spatial_ts1_reduced, exclude_pops = [0])["sample_node_list"]


        ### CHECKS TO MAKE SURE IDs ARE STILL LINING UP ###
        #check to make sure:
        #(1) the focal nodes are coming from the same individual
        #(2) the individual is the correct one (ADMIX_INDIV). This also checks that the individual id hasn't
        #    changed during the simplification step.
        problem_nodes = [node for node in focal_nodes_postsimpl if spatial_ts1_reduced.node(node).individual != ADMIX_INDIV]

        if problem_nodes:
            raise ValueError(
                f"Selected admix nodes {problem_nodes} do not belong to "
                f"focal individual {ADMIX_INDIV}"
            )
        ############

        sample_list_rename = np.unique(spatial_ts1_reduced.nodes_individual[np.where(spatial_ts1_reduced.nodes_flags == 1)])

        indiv_list = []
        node_list = []
        x_list = []
        y_list = []


        for INDIV in sample_list_rename:
            for NODE in spatial_ts1_reduced.individual(INDIV).nodes:
                indiv_list.append(INDIV)
                node_list.append(NODE)
                x_list.append(spatial_ts1_reduced.individual(INDIV).location[0])
                y_list.append(spatial_ts1_reduced.individual(INDIV).location[1])

        indiv_info_dict = {
            "indiv_list":indiv_list,
            "node_list":node_list,
            "x_list":x_list,
            "y_list":y_list
        }


        ts_dict = {
            "indiv": ADMIX_INDIV,
            "tree": spatial_ts1_reduced,
            "node_map": node_map,
            "focal_nodes_postsimpl": focal_nodes_postsimpl,
            "focal_nodes_presimpl": reverse_map[focal_nodes_postsimpl],
            "indiv_info_dict":indiv_info_dict
        }

        simplify_tree_list.append(ts_dict)

    column_names = ['node_id_postsimpl', 'indiv_id', 'ped_id', 'p1_ped_id', 'p2_ped_id']
    output_path = Path(args.output_file_path) / (f"{args.output_prefix}_focal_inds_parent_info.txt")
    with open(output_path, 'w') as file:
        # Write the rows, aligning each list as a column
        file.write("\t".join(column_names) + "\n")
        for ITER in simplify_tree_list:
            for NODE in ITER['focal_nodes_postsimpl']:
                focal_indiv = ITER['tree'].individual(ITER['tree'].node(NODE).individual)
                #focal_indiv = ITER['tree'].node(NODE).individual

                #if focal_indiv != ITER['tree'].node(NODE).individual(ITER['tree'].node(NODE).individual).id:
                #    raise ValueError("The 2 extraction approaches are not the same!")

                #print('stored indiv', ITER['indiv'])
                #print('indiv id from tree', ITER['indiv'])
                if focal_indiv.id != ITER['indiv']:
                    raise ValueError("Unequal indiv ids")

                row = f"{NODE}\t{focal_indiv.id}\t{focal_indiv.metadata['pedigree_id']}\t{focal_indiv.metadata['pedigree_p1']}\t{focal_indiv.metadata['pedigree_p2']}\n"
                file.write(row)


    #Locations of individuals:
    #this file outputs the locations for all individuals/nodes for each round of processing for a focal individual. This is to ensure that
    #the node ids line up with the trees for each focal individual's processing. Note that the focal individual's individual is is included
    #in this file (i.e., it is both focal_indiv_id and indiv_id). I'm not 100% what this location means (I should check on this); it is some
    #residual location passed along from when its ancestors where moved from the main population to the pedigreed population. Currently, I
    #would ignore this.

    #column descriptions:
    #focal_indiv_id: the focal non-georeferenced individual
    #indiv_id: id of individual
    #node_id_postsimpl: node id (after simplification) for the individual specified in indiv_id
    column_names = ['focal_indiv_id', 'indiv_id', 'node_id_postsimpl', 'x', 'y']

    output_path = Path(args.output_file_path) / (f"{args.output_prefix}_location_info.txt")
    with open(output_path, 'w') as file:
        # Write the rows, aligning each list as a column
        file.write("\t".join(column_names) + "\n")

        for INDIV_DICT in simplify_tree_list:
            for i,x in enumerate(INDIV_DICT['indiv_info_dict']['indiv_list']):
                row = f"{INDIV_DICT['indiv']}\t{INDIV_DICT['indiv_info_dict']['indiv_list'][i]}\t{INDIV_DICT['indiv_info_dict']['node_list'][i]}\t{INDIV_DICT['indiv_info_dict']['x_list'][i]}\t{INDIV_DICT['indiv_info_dict']['y_list'][i]}\n"
                file.write(row)

    ### OUTPUTTING TREES AS NEWICKS ###
    if args.tree_output in ("separate", "both"):
        for ITER in simplify_tree_list:
            for TREE in range(0, ITER['tree'].num_trees, args.tree_sep):
                output_tree_path = Path(args.output_tree_path) / (f"{args.output_prefix}_INDIV{str(ITER['indiv'])}_TREE{str(TREE)}.nwk")
                with open(output_tree_path, "w") as f:
                    f.write(ITER['tree'].at_index(TREE).as_newick())

    if args.tree_output in ("gzipped", "both"):
        for ITER in simplify_tree_list:
            output_tree_path = Path(args.output_tree_path) / (f"{args.output_prefix}_INDIV{str(ITER['indiv'])}.nwk.gz")
            with gzip.open(output_tree_path, "wt", encoding="utf-8") as f:
                for TREE in range(0, ITER['tree'].num_trees, args.tree_sep):
                    f.write(ITER['tree'].at_index(TREE).as_newick() + "\n")



    #OUTPUTTING TREE SPAN AND COORDINATE INFO
    column_names = ['focal_indiv', 'tree_index', 'left_interval', 'right_interval', 'span']
    output_path = Path(args.output_file_path) / (f"{args.output_prefix}_tree_interval_info.txt")
    with open(output_path, 'w') as file:
        file.write("\t".join(column_names) + "\n")
        for ITER in simplify_tree_list:
            for TREE in range(0, ITER['tree'].num_trees, args.tree_sep):
                left_interv = ITER['tree'].at_index(TREE).interval.left
                right_interv = ITER['tree'].at_index(TREE).interval.right
                span = ITER['tree'].at_index(TREE).span

                row = f"{ITER['indiv']}\t{TREE}\t{left_interv}\t{right_interv}\t{span}\n"
                file.write(row)


    column_names = ['focal_node_id_presimpl', 'ancestor', 'left_interval', 'right_interval']
    output_path = Path(args.output_file_path) / (f"{args.output_prefix}_ancestry_intervals.txt")
    with open(output_path, 'w') as file:
        file.write("\t".join(column_names) + "\n")
        for index in range(len(anc_node_info_list)):
            focal_node = anc_node_info_list[index]['focal_node']

            for index1 in range(len(anc_node_info_list[index]['ancestor_list'])):
                ancestor = anc_node_info_list[index]['ancestor_list'][index1]
                left_interv = anc_node_info_list[index]['left_interv_list'][index1]
                right_interv = anc_node_info_list[index]['right_interv_list'][index1]

                row = f"{focal_node}\t{ancestor}\t{left_interv}\t{right_interv}\n"
                #row = f"{node_map[focal_node]}\t{focal_node}\t{ancestor}\t{left_interv}\t{right_interv}\n"
                file.write(row)

    # Output file that contains individual id with indiv id and new and old ids
    column_names = ["indiv_id", "node_id_postsimpl", "node_id_presimpl"]
    output_path = Path(args.output_file_path) / (f"{args.output_prefix}_focal_inds_node_map.txt")
    with open(output_path, "w") as file:
        file.write("\t".join(column_names) + "\n")

        for ITER in simplify_tree_list:
            new_nodes = ITER["focal_nodes_postsimpl"]
            old_nodes = ITER["focal_nodes_presimpl"]

            if len(new_nodes) != len(old_nodes):
                raise ValueError(f"Node mapping lengths differ for individual {ITER['indiv']}")

            for new_node, old_node in zip(new_nodes, old_nodes):
                row = f"{ITER['indiv']}\t{new_node}\t{old_node}\n"
                file.write(row)

if __name__ == "__main__":
    args = parser.parse_args()
    main(args)



##################################
### CODE NOT CURRRENTLY IN USE ###
##################################

#(THIS WILL EVENTUALLY BE DELETED ...)

'''
column_names = ['focal_node_simplified', 'focal_node_original', 'ancestor', 'left_interval', 'right_interval']

with open(args.output_file_path + args.output_prefix + '_' + 'node_founder_ancestry_intervals_missingindiv' + '.txt', 'w') as file:
    file.write("\t".join(column_names) + "\n")
    for index in range(len(anc_node_info_list)):
        focal_node = anc_node_info_list[index]['focal_node']

        for index1 in range(len(anc_node_info_list[index]['ancestor_list'])):
            ancestor = anc_node_info_list[index]['ancestor_list'][index1]
            left_interv = anc_node_info_list[index]['left_interv_list'][index1]
            right_interv = anc_node_info_list[index]['right_interv_list'][index1]

            row = f"{node_map[focal_node]}\t{focal_node}\t{ancestor}\t{left_interv}\t{right_interv}\n"
            file.write(row)


    anc_node_info_list = []
    for ADMIX_NODE in admix_nodes_original:
        anc_list = extract_ped_founder_anc(ts = spatial_ts1, focal_node = ADMIX_NODE)
        anc_node_info_list.append(ancestry_intervals(ts = spatial_ts1, anc_list = anc_list, focal_node = ADMIX_NODE))




    #individuals in original population
    non_admixed_indiv_array = sample_indiv_array[~np.isin(sample_indiv_array, admix_indivs)]


    # 2B: individuals in original population
    subsample_sample = np.random.choice(non_admixed_indiv_array, size=args.sample_size, replace=False)
    subsamp_node_list = [spatial_ts1.individual(x).nodes for x in subsample_sample]


    #2C: simplify (and request node map)
    spatial_ts1_reduced,node_map = spatial_ts1.simplify(np.concatenate([np.concatenate(subsamp_node_list), admix_nodes_original]),
                                            map_nodes = True)


    #admix_nodes = spatial_ts2_reduced[0].samples(population = 0)
    admix_nodes = extract_samples_multipop(spatial_ts1_reduced, exclude_pops = [0])["sample_node_list"]

    column_names = ['node', 'indiv', 'ped_id', 'p1_ped_id', 'p2_ped_id']

    with open(args.output_file_path + args.output_prefix + '_admix_node_info.txt', 'w') as file:
        # Write the rows, aligning each list as a column
        file.write("\t".join(column_names) + "\n")
        for NODE in admix_nodes:
            focal_indiv = spatial_ts1_reduced.individual(spatial_ts1_reduced.node(NODE).individual)
            row = f"{NODE}\t{focal_indiv.id}\t{focal_indiv.metadata['pedigree_id']}\t{focal_indiv.metadata['pedigree_p1']}\t{focal_indiv.metadata['pedigree_p2']}\n"
            file.write(row)


    ### EXTRACTING INFO ABOUT INDIVIDUALS (COORDINATES AND NODES ASSOCIATED WITH EACH INDIV) ###
    sample_list_rename = np.unique(spatial_ts1_reduced.nodes_individual[np.where(spatial_ts1_reduced.nodes_flags == 1)])

    indiv_list = []
    node_list = []
    x_list = []
    y_list = []

    for INDIV in sample_list_rename:
        for NODE in spatial_ts1_reduced.individual(INDIV).nodes:
            indiv_list.append(INDIV)
            node_list.append(NODE)
            x_list.append(spatial_ts1_reduced.individual(INDIV).location[0])
            y_list.append(spatial_ts1_reduced.individual(INDIV).location[1])

    column_names = ['indiv', 'node', 'x', 'y']

    with open(args.output_file_path + args.output_prefix + '_indiv_node_info.txt', 'w') as file:
        # Write the rows, aligning each list as a column
        file.write("\t".join(column_names) + "\n")
        for i,x in enumerate(indiv_list):
                row = f"{indiv_list[i]}\t{node_list[i]}\t{x_list[i]}\t{y_list[i]}\n"
                file.write(row)


    ### OUTPUTTING TREES AS NEWICKS ###
    if args.tree_output in ("separate", "both"):
        for TREE in range(0, spatial_ts1_reduced.num_trees, args.tree_sep):
            with open(args.output_tree_path + args.output_prefix + str(TREE) + ".nwk", "w") as f:
                f.write(spatial_ts1_reduced.at_index(TREE).as_newick())

    if args.tree_output in ("gzipped", "both"):
        with gzip.open(args.output_tree_path + args.output_prefix + ".nwk.gz", "wt", encoding="utf-8") as f:
            for TREE in range(0, spatial_ts1_reduced.num_trees, args.tree_sep):
                f.write(spatial_ts1_reduced.at_index(TREE).as_newick() + "\n")
            #print('exported' + str(TREE))



    #OUTPUTTING TREE SPAN AND COORDINATE INFO
    column_names = ['tree_index', 'left_interval', 'right_interval', 'span']
    with open(args.output_file_path + 'tree_interval_info.txt', 'w') as file:
        file.write("\t".join(column_names) + "\n")
        for TREE in range(0, spatial_ts1_reduced.num_trees, args.tree_sep):
            left_interv = spatial_ts1_reduced.at_index(TREE).interval.left
            right_interv = spatial_ts1_reduced.at_index(TREE).interval.right
            span = spatial_ts1_reduced.at_index(TREE).span

            row = f"{TREE}\t{left_interv}\t{right_interv}\t{span}\n"
            file.write(row)


    column_names = ['focal_node_simplified', 'focal_node_original', 'ancestor', 'left_interval', 'right_interval']

    with open(args.output_file_path + 'node_founder_ancestry_intervals.txt', 'w') as file:
        file.write("\t".join(column_names) + "\n")
        for index in range(len(anc_node_info_list)):
            focal_node = anc_node_info_list[index]['focal_node']

            for index1 in range(len(anc_node_info_list[index]['ancestor_list'])):
                ancestor = anc_node_info_list[index]['ancestor_list'][index1]
                left_interv = anc_node_info_list[index]['left_interv_list'][index1]
                right_interv = anc_node_info_list[index]['right_interv_list'][index1]

                row = f"{node_map[focal_node]}\t{focal_node}\t{ancestor}\t{left_interv}\t{right_interv}\n"
                file.write(row)



    #admixed individuals
    admix_indivs = list(set([spatial_ts1.node(i).individual for i in admix_nodes_original])) #individuals in the admixed pop
'''


'''
    column_names = ['focal_indiv', 'focal_node_simplified', 'focal_node_original', 'ancestor', 'left_interval', 'right_interval']

    with open(args.output_file_path + args.output_prefix + '_' + 'node_founder_ancestry_intervals_missingindiv' + '.txt', 'w') as file:
        file.write("\t".join(column_names) + "\n")
        for index in range(len(anc_node_info_list)):
            focal_node = anc_node_info_list[index]['focal_node']

            for index1 in range(len(anc_node_info_list[index]['ancestor_list'])):
                ancestor = anc_node_info_list[index]['ancestor_list'][index1]
                left_interv = anc_node_info_list[index]['left_interv_list'][index1]
                right_interv = anc_node_info_list[index]['right_interv_list'][index1]

                for ITER in simplify_tree_list:
                    row = f"{ITER['indiv']}\t{ITER['node_map'][focal_node]}\t{focal_node}\t{ancestor}\t{left_interv}\t{right_interv}\n"
                    file.write(row)


column_names = ['node', 'indiv', 'ped_id', 'p1_ped_id', 'p2_ped_id']

with open(args.output_file_path + args.output_prefix + '_admix_node_info.txt', 'w') as file:
     # Write the rows, aligning each list as a column
    file.write("\t".join(column_names) + "\n")
    for ITER in simplify_tree_list:
        for NODE in ITER['admix_nodes']:
            focal_indiv = ITER['tree'].individual(ITER['tree'].node(NODE).individual)
            if focal_indiv.id != ITER['indiv']:
                raise ValueError("Unequal indiv ids")

            row = f"{NODE}\t{focal_indiv.id}\t{focal_indiv.metadata['pedigree_id']}\t{focal_indiv.metadata['pedigree_p1']}\t{focal_indiv.metadata['pedigree_p2']}\n"
            file.write(row)


sample_list_rename = np.unique(spatial_ts1_reduced.nodes_individual[np.where(spatial_ts1_reduced.nodes_flags == 1)])



column_names = ['focal_indiv', 'indiv', 'node', 'x', 'y']

with open(args.output_file_path + args.output_prefix + '_indiv_node_info.txt', 'w') as file:
    # Write the rows, aligning each list as a column
    file.write("\t".join(column_names) + "\n")

    for INDIV_DICT in simplify_tree_list:
        for i,x in enumerate(INDIV_DICT['indiv_info_dict']['indiv_list']):
            row = f"{INDIV_DICT['indiv']}{INDIV_DICT['indiv_info_dict']['indiv_list'][i]}\t{INDIV_DICT['indiv_info_dict']['node_list'][i]}\t{INDIV_DICT['indiv_info_dict']['x_list'][i]}\t{INDIV_DICT['indiv_info_dict']['y_list'][i]}\n"
            file.write(row)




### OUTPUTTING TREES AS NEWICKS ###
if args.tree_output in ("separate", "both"):
    for TREE in range(0, spatial_ts1_reduced.num_trees, args.tree_sep):
        with open(args.output_tree_path + args.output_prefix + str(TREE) + ".nwk", "w") as f:
            f.write(spatial_ts1_reduced.at_index(TREE).as_newick())

if args.tree_output in ("gzipped", "both"):
    with gzip.open(args.output_tree_path + args.output_prefix + ".nwk.gz", "wt", encoding="utf-8") as f:
        for TREE in range(0, spatial_ts1_reduced.num_trees, args.tree_sep):
            f.write(spatial_ts1_reduced.at_index(TREE).as_newick() + "\n")






### EXTRACTING INFO ABOUT INDIVIDUALS (COORDINATES AND NODES ASSOCIATED WITH EACH INDIV) ###
    sample_list_rename = np.unique(spatial_ts1_reduced.nodes_individual[np.where(spatial_ts1_reduced.nodes_flags == 1)])

    indiv_list = []
    node_list = []
    x_list = []
    y_list = []

    for INDIV in sample_list_rename:
        for NODE in spatial_ts1_reduced.individual(INDIV).nodes:
            indiv_list.append(INDIV)
            node_list.append(NODE)
            x_list.append(spatial_ts1_reduced.individual(INDIV).location[0])
            y_list.append(spatial_ts1_reduced.individual(INDIV).location[1])

    column_names = ['indiv', 'node', 'x', 'y']

    with open(args.output_file_path + args.output_prefix + '_indiv_node_info.txt', 'w') as file:
        # Write the rows, aligning each list as a column
        file.write("\t".join(column_names) + "\n")
        for i,x in enumerate(indiv_list):
                row = f"{indiv_list[i]}\t{node_list[i]}\t{x_list[i]}\t{y_list[i]}\n"
                file.write(row)



for ADMIX_INDIV in admix_indivs:

    focal_nodes = spatial_ts1.individuals(ADMIX_INDIV).nodes

    #2C: simplify (and request node map)
    spatial_ts1_reduced,node_map = spatial_ts1.simplify(np.concatenate([np.concatenate(subsamp_node_list), focal_nodes]),
                                            map_nodes = True)


    #admix_nodes = spatial_ts2_reduced[0].samples(population = 0)
    admix_nodes = extract_samples_multipop(spatial_ts1_reduced, exclude_pops = [0])["sample_node_list"]

    column_names = ['focal_indiv', 'node', 'indiv', 'ped_id', 'p1_ped_id', 'p2_ped_id']

    with open(args.output_file_path + args.output_prefix + '_admix_node_info.txt', 'w') as file:
        # Write the rows, aligning each list as a column
        file.write("\t".join(column_names) + "\n")
        for NODE in admix_nodes:
            focal_indiv = spatial_ts1_reduced.individual(spatial_ts1_reduced.node(NODE).individual)
            row = f"{NODE}\t{focal_indiv.id}\t{focal_indiv.metadata['pedigree_id']}\t{focal_indiv.metadata['pedigree_p1']}\t{focal_indiv.metadata['pedigree_p2']}\n"
            file.write(row)


    ### EXTRACTING INFO ABOUT INDIVIDUALS (COORDINATES AND NODES ASSOCIATED WITH EACH INDIV) ###
    sample_list_rename = np.unique(spatial_ts1_reduced.nodes_individual[np.where(spatial_ts1_reduced.nodes_flags == 1)])

    indiv_list = []
    node_list = []
    x_list = []
    y_list = []

    for INDIV in sample_list_rename:
        for NODE in spatial_ts1_reduced.individual(INDIV).nodes:
            indiv_list.append(INDIV)
            node_list.append(NODE)
            x_list.append(spatial_ts1_reduced.individual(INDIV).location[0])
            y_list.append(spatial_ts1_reduced.individual(INDIV).location[1])

    column_names = ['indiv', 'node', 'x', 'y']

    with open(args.output_file_path + args.output_prefix + '_indiv_node_info.txt', 'w') as file:
        # Write the rows, aligning each list as a column
        file.write("\t".join(column_names) + "\n")
        for i,x in enumerate(indiv_list):
                row = f"{indiv_list[i]}\t{node_list[i]}\t{x_list[i]}\t{y_list[i]}\n"
                file.write(row)


    ### OUTPUTTING TREES AS NEWICKS ###
    if args.tree_output in ("separate", "both"):
        for TREE in range(0, spatial_ts1_reduced.num_trees, args.tree_sep):
            with open(args.output_tree_path + args.output_prefix + str(TREE) + ".nwk", "w") as f:
                f.write(spatial_ts1_reduced.at_index(TREE).as_newick())

    if args.tree_output in ("gzipped", "both"):
        with gzip.open(args.output_tree_path + args.output_prefix + ".nwk.gz", "wt", encoding="utf-8") as f:
            for TREE in range(0, spatial_ts1_reduced.num_trees, args.tree_sep):
                f.write(spatial_ts1_reduced.at_index(TREE).as_newick() + "\n")
            #print('exported' + str(TREE))



    #OUTPUTTING TREE SPAN AND COORDINATE INFO
    column_names = ['tree_index', 'left_interval', 'right_interval', 'span']
    with open(args.output_file_path + 'tree_interval_info.txt', 'w') as file:
        file.write("\t".join(column_names) + "\n")
        for TREE in range(0, spatial_ts1_reduced.num_trees, args.tree_sep):
            left_interv = spatial_ts1_reduced.at_index(TREE).interval.left
            right_interv = spatial_ts1_reduced.at_index(TREE).interval.right
            span = spatial_ts1_reduced.at_index(TREE).span

            row = f"{TREE}\t{left_interv}\t{right_interv}\t{span}\n"
            file.write(row)


    column_names = ['focal_node_simplified', 'focal_node_original', 'ancestor', 'left_interval', 'right_interval']

    with open(args.output_file_path + 'node_founder_ancestry_intervals_missingindiv' + ADMIX_INDIV + '.txt', 'w') as file:
        file.write("\t".join(column_names) + "\n")
        for index in range(len(anc_node_info_list)):
            focal_node = anc_node_info_list[index]['focal_node']

            for index1 in range(len(anc_node_info_list[index]['ancestor_list'])):
                ancestor = anc_node_info_list[index]['ancestor_list'][index1]
                left_interv = anc_node_info_list[index]['left_interv_list'][index1]
                right_interv = anc_node_info_list[index]['right_interv_list'][index1]

                row = f"{node_map[focal_node]}\t{focal_node}\t{ancestor}\t{left_interv}\t{right_interv}\n"
                file.write(row)
'''







########################################
### CODE THAT'S NOT CURRENTLY IN USE ###
########################################

### EXTRACTING INFO ABOUT THE TREES (COORDINATES OF ROOT AND GENOME SPAN) ###
'''
tree_index_list = []
tree_position_list = []
node_list = []
x_list = []
y_list = []
span_list = []

for i,tree_pos in enumerate(range(0, spatial_ts2_reduced[0].num_trees, args.tree_sep)):
    tree_index_list.append(i)
    tree_position_list.append(tree_pos)
    x_list.append(spatial_ts2_reduced[0].individual(spatial_ts2_reduced[0].node(spatial_ts2_reduced[0].at_index(tree_pos).root).individual).location[0])
    y_list.append(spatial_ts2_reduced[0].individual(spatial_ts2_reduced[0].node(spatial_ts2_reduced[0].at_index(tree_pos).root).individual).location[1])
    node_list.append(spatial_ts2_reduced[0].node(spatial_ts2_reduced[0].at_index(tree_pos).root).id)
    span_list.append(spatial_ts2_reduced[0].at_index(tree_pos).span)

column_names = ['tree_index', 'tree_position', 'node', 'x', 'y', 'span']


with open(args.output_file_path + args.output_prefix + '_tree_position_info.txt', 'w') as file:
    # Write the rows, aligning each list as a column
    file.write("\t".join(column_names) + "\n")
    for i,x in enumerate(tree_index_list):
            row = f"{tree_index_list[i]}\t{tree_position_list[i]}\t{node_list[i]}\t{x_list[i]}\t{y_list[i]}\t{span_list[i]}\n"
            file.write(row)
'''




'''
def ancestor_at_time_multitree_nodes(ts, nodes, time):
    tree_index_list = []
    node_list = []
    sample_node_list = []

    for NODE in nodes:
        output = ancestor_at_time_multitree(ts = ts, node = NODE, time = time)
        tree_index_list.extend(output['tree_index'])
        node_list.extend(output['node'])
        sample_node_list.extend(ts.num_trees * [NODE])

    return({
        "tree_index": tree_index_list,
        "ancestor_node": node_list,
        "sample_node": sample_node_list
    })


#COMMENT OUT?
def ancestor_at_time_multitree(ts, node, time):
    node_list = []
    for TREE in ts.trees():
        node_list.append(ancestor_at_time(tree = TREE, node = node, time = time))

    #[i[0] for i in node_list]
    return({"tree_index": [i for i in range(0, ts.num_trees, 1)],
            "node":node_list})

#COMMENT OUT?
def ancestor_at_time(tree, node, time):
    node_collector = [node]
    node_age = [tree.time(node)]
    new_time = tree.time(node)
    while new_time < time:
        new_node = tree.parent(node_collector[-1])
        new_time = tree.time(new_node)

        if new_time <= time:
            node_collector.append(new_node)
            node_age.append(new_time)

        #return(
    #   [node_collector[i] if abs(age - time) <= 1e-15 else -9999 for i, age in enumerate(node_age)]
    #)
    if abs(node_age[-1] - time) <= 1e-15:
        return(node_collector[-1])
    else:
        return(-9999)

#COMMENT OUT
def combine_elements_index(val):
    new_val = []
    ind_list = []
    index = 0
    for i,VAL in enumerate(val):
        if index == 0:
            left_ind = i
            new_val.append(VAL)
            index += 1

        if i == (len(val) - 1):
            ind_list.append([left_ind, i])
        elif VAL != val[i + 1]:
            ind_list.append([left_ind, i])
            index = 0
    return(
        {"val":new_val,
         "ind_list":ind_list}
    )
'''


'''
### STEP 3: FOR EACH ADMIXED INDIVIDUAL, IDENTIFY ITS ANCESTRY (SPECIFIC NODE AND INDIV)
###         FROM THE MAIN POPULATION
###        - THIS IS IMPORTANT BECAUSE IT WILL TELL US FROM WHERE IN THE ORIGINAL POPULATION
###          EACH ADMIXED ANCESTORS IS FROM
# 3A: IDENTIFY ADMIXED NODES
# 3B: IDENTIFY THE TIMING OF THE INDIVIDUALS
# 3C: FIND THE ANCESTOR NODE AND IN EACH TREE FOR EACH SAMPLE NODE
# 3D: CONNECT EACH ANCESTOR NODE TO THE PEDIGREE ID


# 3A
admix_nodes_simp = spatial_ts1_reduced.samples(population=1)

# 3B
parent_pop2_gen = max(spatial_ts1_reduced.nodes_time[spatial_ts1_reduced.nodes_population == 1]) + 1.0

# 3C
# outputs a dictionary containing 3 lists:
# - sample node
# - ancestor node
# - tree index
ancestor_node_info = ancestor_at_time_multitree_nodes(ts = spatial_ts1_reduced,
                                        nodes = admix_nodes_simp,
                                        time = parent_pop2_gen)

# 3D
parent_pedigree_id_list = []
for NODE in ancestor_node_info['ancestor_node']:
    parent_pedigree_id_list.append(spatial_ts1_reduced.individual(spatial_ts1_reduced.node(NODE).individual).metadata['pedigree_id'])
'''


### STEP 4: SIMPLIFY AGAIN, BUT REMOVE UNARY NODES
# 4A: SIMPLIFY
# 4B: CHECK NODE MAPPING

#notes: for some reason, this results in the reduction of tree count while keep_unary does not. I'm not sure
#why the keep unary doesn't reduce the number of trees.

#This is the final set of trees that are output, so we need to know the ancestry of each admixed individual in
#this tree sequence even though though we have simplified out the relevant unary nodes. I'm currently sure how
#to perform this simplification, but retain the relevant unary nodes from the p1 population ... The way that I
#doing this seems much complicated, convoluted, and slow ...

#DELETE
#spatial_ts2_reduced = spatial_ts1_reduced.simplify(map_nodes = True)


'''
#because simplification can result in the changing of node and indiv IDs, check to see that the ids of the
#admixed nodes have not changed. I am doing this by setting the map_nodes to True, and then checking that
#the admixed node IDs from the previous tree sequence (prior to the most recent simplification) are equal
#to the value of the entry in the map_nodes array at the index value equal to the previous node id. I'm
#not 100% that this is correct ... but I think it is?
node_id_switch_check = (admix_nodes_simp==spatial_ts2_reduced[1][admix_nodes_simp]).all()
if not node_id_switch_check:
    raise TypeError("Not all nodes are the same.")
else:
    print('no changes')
'''

### STEP 5: COMBINE INFORMATION ACROSS TREES TO FIGURE OUT FOR A SPECIFIC SAMPLE NODE, WHAT
###         COORDINATES IN THE GENOME DOES IT HAVE ANCESTRY FROM A PARTICULAR INDIVIDUAL
#- this is one (very convoluted way) of mapping the ancestry info from one set of trees to
#  another set of (more simplified) trees
#- an easier way would be to just directly get this info from edge information? Like instead
#  of going tree-by-tree in this and previous steps, just getting the coordinates from the
#  edge table and not bothering with the tree-wise extraction and combination

#THE REASON THIS IS SO SLOW IS THE at_index CALL IS SEEMINGLY VERY VERY VERY SLOW

'''
node_array = np.array(ancestor_node_info['sample_node'])
sample_nodes = np.unique(node_array)

sample_node_list = []
ancestor_node = []
ped_id = []
ancestor_ped_id = []
left_span = []
right_span = []


for NODE in sample_nodes:
    print('started ' + str(NODE))
    indices = np.where(node_array == NODE)[0]
    ancestor_node_array = np.array(ancestor_node_info['ancestor_node'])[indices]

    test_combine = combine_elements_index(ancestor_node_array)

    for i,ANC_NODE in enumerate(test_combine['val']):
        sample_node_list.append(NODE)
        ancestor_node.append(ANC_NODE)
        ped_id.append(spatial_ts1_reduced.individual(spatial_ts1_reduced.node(NODE).individual).metadata['pedigree_id'])
        ancestor_ped_id.append(spatial_ts1_reduced.individual(spatial_ts1_reduced.node(ANC_NODE).individual).metadata['pedigree_id'])

        left_int_info = spatial_ts1_reduced.at_index(test_combine['ind_list'][i][0]).interval
        left_span.append(left_int_info[0])

        right_int_info = spatial_ts1_reduced.at_index(test_combine['ind_list'][i][1]).interval
        right_span.append(right_int_info[1])
        print('finished ' + str(i))
'''



'''
### STEP 7: FOR EACH ADMIXED NODE, IDENTIFY THE ANCESTRAL NODE AND ASSOCIATED PEDIGREE
###         ID IN EACH TREE.
# - this is necessary because the tree simplified with unary nodes retained (which is
#   necessary for extracting relevant ancestry info) has more local trees than the
#   ts simplified without retaining unary nodes. The task is using the information
#   from the unary nodes trees regarding ancestry to figure out the ancestry in the
#   ts where simplification involved removing unary nodes

final_sample_node_list = []
final_ancestor_node_list = []
final_ped_id_list = []
tree_index_list = []

sample_node_array = np.unique(np.array(sample_node_list))


for SAMPLE_NODE in sample_node_array:

    indices = np.where(np.array(sample_node_list) == SAMPLE_NODE)[0] # indices associated with sample node
    left_span_sub = np.array(left_span)[indices] #the left span for each tree for sample node
    right_span_sub = np.array(right_span)[indices] #theright span for each for sample node
    ancestor_node_sub = np.array(ancestor_node)[indices] #the ancestor node for sample node
    ped_id_sub = np.array(ancestor_ped_id)[indices]

    for index, (LEFT, RIGHT) in enumerate(zip(new_simp_left_span, new_simp_right_span)):
        #left_ind = np.max(np.where(left_span_sub <= new_simp_left_span[25000]))
        #right_ind = np.min(np.where(right_span_sub >= new_simp_right_span[25000]))
        left_ind = np.max(np.where(left_span_sub <= LEFT))
        right_ind = np.min(np.where(right_span_sub >= RIGHT))

        if left_ind != right_ind:
            print('left_ind = ' + str(left_ind))
            print('right_ind = ' + str(right_ind))

            final_sample_node_list.append(SAMPLE_NODE)
            final_ancestor_node_list.append(-9999)

            if ped_id_sub[left_ind] == ped_id_sub[right_ind]:
                final_ped_id_list.append(ped_id_sub[left_ind])
            else:
                final_ped_id_list.append(-9999)
        else:
            final_sample_node_list.append(SAMPLE_NODE)
            final_ancestor_node_list.append(ancestor_node_sub[left_ind])
            final_ped_id_list.append(ped_id_sub[left_ind])

            #raise ValueError("Indices do not match")
        tree_index_list.append(index)


    print('finished ' + str(SAMPLE_NODE))
'''


'''
column_names = ['tree_index', 'sample_node', 'ancestor_node', 'ancestor_ped_id']

with open(output_file_path + 'ancester_info_simplified_tree.txt', 'w') as file:
    # Write the rows, aligning each list as a column
    file.write("\t".join(column_names) + "\n")
    for i,x in enumerate(tree_index_list):
            row = f"{tree_index_list[i]}\t{final_sample_node_list[i]}\t{final_ancestor_node_list[i]}\t{final_ped_id_list[i]}\n"
            file.write(row)

'''