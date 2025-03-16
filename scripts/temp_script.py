import dendropy
from pandas import Series, DataFrame
import numpy as np


'''
AnnArtiges
Temporary script to retrieve information in newick files.
For now we retrieve the distances between each leaf in a tree and the paralogs and common leaves in two trees.
'''


def make_distance_matrix(distances, sort=False):
    distances = distances.reset_index().pivot_table(index=0, columns=1, values='distance')
    distances.index.name=None
    distances.columns.name=None
    distances = distances.fillna(distances.T)

    missing_rows = [i for i in distances.index if i not in distances.columns]
    for row in missing_rows:
        distances[row] = distances.loc[row]

    missing_columns = [c for c in distances.columns if c not in distances.index]
    for column in missing_columns:
        distances.loc[column] = distances[column]
    
    if sort:
        distances = distances.sort_index().sort_index(axis=1)
    return distances.fillna(0)



def d(start, end, distances):
    return distances.loc[(start, end), 'distance']



def get_distance_dataframe(tree):
    pdm = tree.phylogenetic_distance_matrix()
    distances = [[*name, distance] for name, distance in zip(pdm.distinct_taxon_pair_iter(), pdm.distances())]
    distances = DataFrame(distances)
    
    def get_label(node):
        return node.label
    
    distances[0] = distances[0].map(get_label)
    distances[1] = distances[1].map(get_label)
    
    distances = distances.set_index([0, 1])
    distances = distances.rename(columns = {2:'distance'})

    return distances



def get_ratio(distances1, distances2):
    r1 = distances1.loc[common_leaves, common_leaves].sum().sum() / 2
    r2 = distances2.loc[common_leaves, common_leaves].sum().sum() / 2
    return r1/r2

def noisify_distances(distances):
    #Function to use with generated trees => keep the same topology/organisation but with different branch lengths

    noise = np.triu(1 + np.random.random(distances.shape) / 10).round(2)
    noise += noise.T
    return distances * noise


#def temp_paralog(taxalist):  


def get_paralogs(taxa_list1, taxa_list2):
    # Temporary function retrieveing the paralogs and common leaves in two trees.
    # Updates/modifications to implement :  - add the tree corresponding to each paralogs
    #                                       - give new names to paralogs and store original ID somewhere
    #

    common_leaves = []
    temp_leaves = []
    paralogs_species = []
    dico_identifier = {}
    
    CommonLeaves = {}
    Paralogs = {}
    
    for c in taxa_list1:
        taxon = c.split(' ')[0]
        gene = c.split(' ')[1]
        
        if taxon in dico_identifier.keys() :
            dico_identifier[taxon].append(gene)
        else : dico_identifier[taxon] = [gene]

        if taxon in temp_leaves :
            paralogs_species.append(taxon)
            temp_leaves.remove(taxon)
            
        else :
            temp_leaves.append(taxon)

    
    for d in taxa_list2:
        taxon = d.split(' ')[0]
        gene = d.split(' ')[1]
        
        if taxon in dico_identifier.keys() :
            dico_identifier[taxon].append(gene)
        else : dico_identifier[taxon] = [gene]
        
        if taxon in temp_leaves :
            common_leaves.append(taxon)
            temp_leaves.remove(taxon)
            
        elif taxon in common_leaves :
            paralogs_species.append(taxon)
            common_leaves.remove(taxon)
    
    for specie in paralogs_species : Paralogs[specie] = dico_identifier[specie]
    for common in common_leaves : CommonLeaves[common] = dico_identifier[common]
    
    #print(Paralogs['n212'])
    #print(CommonLeaves['n1'])
    
    return paralogs_species, common_leaves




tree = dendropy.Tree.get(path='data/Sample_Gf2984.tre', schema='newick') # or whatever relevant format if not newick

'''list_name = []
for name in tree.phylogenetic_distance_matrix().distinct_taxon_pair_iter():
    if name[0] not in list_name :
        list_name.append(name[0])
    elif name[1] not in list_name :
        list_name.append(name[1])'''


distances1 = get_distance_dataframe(tree)

test = make_distance_matrix(distances1, True)
taxa = test.columns

common_leaves = {}
paralogs = {}
paralogs, common_leaves = get_paralogs(taxa, taxa)