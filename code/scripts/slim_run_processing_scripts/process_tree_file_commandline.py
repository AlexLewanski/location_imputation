########################
### LOADING PACKAGES ###
########################
#import sys
import argparse
import tskit
#import msprime
import numpy as np
import pyslim
#import networkx as nx
#from IPython.display import SVG
#import math
import random

#tree_separation = 50
#ts_upload_path = '/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output/admix_sim/'
#output_file_path = '/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output/admix_sim/'
#output_tree_path = '/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output/admix_sim/' #tree_dir/'


parser = argparse.ArgumentParser(description="Process tree sequence.")

# Positional arguments
parser.add_argument("--tree_sep", "-t", type=int, help="separation for tree sampling")
parser.add_argument("--sample_size", "-s", type=int, help="separation for tree sampling")
parser.add_argument("--ts_upload_path", "-u", help="upload path for processed tree sequence")
parser.add_argument("--ts_name", "-m", help="name of tree sequence to upload")
parser.add_argument("--output_file_path", "-o", help="output path for files")
parser.add_argument("--output_tree_path", "-p", help="output path for newick trees")
parser.add_argument("--output_prefix", "-f", help="output name for newick trees")
parser.add_argument("--recomb_rate", "-r", type=float, help="recombination rate for msprime simulation")
parser.add_argument("--ne", "-n", type=float, help="effective size for msprime simulation")
parser.add_argument("--seed", "-d", type=float, help="random seed")

args = parser.parse_args()


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
    

def ancestor_at_time_multitree(ts, node, time):
    node_list = []
    for TREE in ts.trees():
        node_list.append(ancestor_at_time(tree = TREE, node = node, time = time))
    
    #[i[0] for i in node_list]
    return({"tree_index": [i for i in range(0, ts.num_trees, 1)],
            "node":node_list})

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

def extract_samples_multipop(ts, exclude_pops):
    sample_node_list=[]
    sample_node_pop=[]
    for i in ts.samples():
        if ts.node(i).population not in exclude_pops:
            sample_node_list.append(i)
            sample_node_pop.append(ts.node(i).population)

    return({"sample_node_list": sample_node_list,
            "pop_list": sample_node_pop})



###############################################
### LOADING AND INITIAL PROCESSING OF TREES ###
###############################################

### PROCESSING VARIABLES ###
#tree_separation = 50

#ts_upload_path = '/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output/admix_sim/'
#output_file_path = '/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output/admix_sim/'
#output_tree_path = '/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output/admix_sim/' #tree_dir/'



### STEP 1: INITIAL STEPS
# 1A: LOAD INITIAL TREE SEQUENCE
# 1B: RECAPITATE TREE SEQUENCE
# ***NO SIMPLIFICATION

#1A
spatial_ts = tskit.load(args.ts_upload_path + args.ts_name)
#spatial_ts = tskit.load("/Users/alexlewanski/Documents/michigan_state/research/location_imputation/simulation_output/postburnin_sims/POSTBURNIN_BURNINREP1_RUNTIME10000_DISP1.80_CC0.15_K3_Kreps15.trees")

#1B
#spatial_ts1 = pyslim.recapitate(
#    spatial_ts,
#    recombination_rate = args.recomb_rate,
#    ancestral_Ne = args.ne,
#    random_seed=71110)

#1B
spatial_ts1 = pyslim.recapitate(
    spatial_ts,
    recombination_rate = args.recomb_rate,
    ancestral_Ne = args.ne,
    random_seed=args.seed)


random.seed(args.seed)

### STEP 2: INITIAL SIMPLIFICATION INCLUDING ADMIXED INDIVS + OTHERS
# 2A: IDENTIFY ADMIXED INDIVS
# 2B: IDENTIFY NON-ADMIXED INDIVS TO BE SAMPLED 
# 2C: SIMPLIFY BASED ON ADMIXED AND SAMPLED INDIVS (KEEP UNARY NODES!)

#sample individuals
sample_indiv_array = np.unique(spatial_ts1.nodes_individual[np.where(spatial_ts1.nodes_flags == 1)]) #sample individuals

#2A
#admix_nodes = spatial_ts1.samples(population=2) #nodes in admixed population
admix_nodes_dict = extract_samples_multipop(ts = spatial_ts1, exclude_pops = [1])
admix_nodes=admix_nodes_dict["sample_node_list"]
admix_indivs = list(set([spatial_ts1.node(i).individual for i in admix_nodes])) #individuals in the admixed pop


#individuals in original population
non_admixed_indiv_array = sample_indiv_array[~np.isin(sample_indiv_array, admix_indivs)]


# 2B: individuals in original population
subsample_sample = np.random.choice(non_admixed_indiv_array, size=args.sample_size, replace=False)
subsamp_node_list = [spatial_ts1.individual(x).nodes for x in subsample_sample]


#2C: simplify
spatial_ts1_reduced = spatial_ts1.simplify(np.concatenate([np.concatenate(subsamp_node_list), admix_nodes]),
                                           keep_unary = True)


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


spatial_ts2_reduced = spatial_ts1_reduced.simplify(map_nodes = True)

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
#admix_nodes = spatial_ts2_reduced[0].samples(population = 0)
admix_nodes = extract_samples_multipop(spatial_ts2_reduced[0], exclude_pops = [0])["sample_node_list"]

column_names = ['node', 'indiv', 'ped_id', 'p1_ped_id', 'p2_ped_id']

with open(args.output_file_path + args.output_prefix + '_admix_node_info.txt', 'w') as file:
    # Write the rows, aligning each list as a column
    file.write("\t".join(column_names) + "\n")
    for NODE in admix_nodes:
        focal_indiv = spatial_ts2_reduced[0].individual(spatial_ts2_reduced[0].node(NODE).individual)
        row = f"{NODE}\t{focal_indiv.id}\t{focal_indiv.metadata['pedigree_id']}\t{focal_indiv.metadata['pedigree_p1']}\t{focal_indiv.metadata['pedigree_p2']}\n"
        file.write(row)



### EXTRACTING INFO ABOUT THE TREES (COORDINATES OF ROOT AND GENOME SPAN) ###
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


### EXTRACTING INFO ABOUT INDIVIDUALS (COORDINATES AND NODES ASSOCIATED WITH EACH INDIV) ###
sample_list_rename = np.unique(spatial_ts2_reduced[0].nodes_individual[np.where(spatial_ts2_reduced[0].nodes_flags == 1)])

indiv_list = []
node_list = []
x_list = []
y_list = []

for INDIV in sample_list_rename:
    for NODE in spatial_ts2_reduced[0].individual(INDIV).nodes:
        indiv_list.append(INDIV)
        node_list.append(NODE)
        x_list.append(spatial_ts2_reduced[0].individual(INDIV).location[0])
        y_list.append(spatial_ts2_reduced[0].individual(INDIV).location[1])

column_names = ['indiv', 'node', 'x', 'y']

with open(args.output_file_path + args.output_prefix + '_indiv_node_info.txt', 'w') as file:
    # Write the rows, aligning each list as a column
    file.write("\t".join(column_names) + "\n")
    for i,x in enumerate(indiv_list):
            row = f"{indiv_list[i]}\t{node_list[i]}\t{x_list[i]}\t{y_list[i]}\n"
            file.write(row)


### OUTPUTTING TREES AS NEWICKS ###
for TREE in range(0, spatial_ts2_reduced[0].num_trees, args.tree_sep):
    with open(args.output_tree_path + args.output_prefix + str(TREE) + ".nwk", "w") as f:
        f.write(spatial_ts2_reduced[0].at_index(TREE).as_newick())
        #print('exported' + str(TREE))
