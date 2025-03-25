#!/usr/bin/python
import dendropy
from pandas import Series, DataFrame
import numpy as np

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
    #print(distances)
    
    def get_label(node):
        return node.label
    
    distances[0] = distances[0].map(get_label)
    distances[1] = distances[1].map(get_label)
    
    distances = distances.set_index([0, 1])
    distances = distances.rename(columns = {2:'distance'})

    return distances


#Same function as the get_distance_dataframe but only looks at the trees' topology.
def get_distance_topology(tree):
    pdm = tree.phylogenetic_distance_matrix()
    distances = [[*name, distance] for name, distance in zip(pdm.distinct_taxon_pair_iter(), pdm.distances())]
    
    for pairs in distances:
        pairs[2] = pdm.path_edge_count(taxon1 = pairs[0], taxon2 = pairs[1])
    
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

#Temp function to create tree with same topology but different branch lengths.
#Can serve as control : trees with same topology should have the same pairings (if the branch lengths have not change drastically)
def noisify_distances(distances):
    noise = np.triu(1 + np.random.random(distances.shape) / 10).round(2)
    noise += noise.T
    return distances * noise


#def temp_paralog(taxalist):
    

#Temp function to identify the paralogs in a tree and the common leaves between two trees.
def get_paralogs(taxa_list1, taxa_list2):
    
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
    
    return Paralogs, CommonLeaves




def get_paralogs_variant(taxa1, taxa2, common_leaves, paralogs):
    
    dict_T1 = {"Specie":[], "Identifier":[], "Alias":[], "Distance":[]}
    dict_T2 = {"Specie":[], "Identifier":[], "Alias":[], "Distance":[]}
    
    for c in taxa1:
        taxon = c.split(' ')[0]
        gene = c.split(' ')[1]
        dict_T1["Specie"].append(taxon)
        dict_T1["Identifier"].append(gene)
        dict_T1["Alias"].append(' ')
        dict_T1["Distance"].append(0)
    
    for d in taxa2:
        taxon = d.split(' ')[0]
        gene = d.split(' ')[1]
        dict_T2["Specie"].append(taxon)
        dict_T2["Identifier"].append(gene)
        dict_T2["Alias"].append(' ')
        dict_T2["Distance"].append(0)
        
    temp_common = np.intersect1d(dict_T1["Specie"],dict_T2["Specie"])
    for element in temp_common:
        if dict_T1["Specie"].count(element) == 1 and dict_T2["Specie"].count(element) == 1:
            common_leaves.append(element)
        else : paralogs.append(element)
    
    return dict_T1, dict_T2

def get_SumDistance(common_leaves, matrix, dic):
    sum_CL = 0
    i = 1
    
    for CL in common_leaves[:-1]:
        leaf1 = CL+" "+dic["Identifier"][dic["Specie"].index(CL)]
        col = matrix[leaf1]
        for cl in common_leaves[i:]:
            leaf2 = cl+" "+dic["Identifier"][dic["Specie"].index(cl)]
            sum_CL += col[leaf2]
        i+=1


'''
#######
#######  NOPE MARCHE PAS, MODIFICATION À APPORTER POUR AVOIR LA DISTANCE POUR CHAQUE PARALOG ET PAS JUSTE LE TOTAL
#######                   ==> ajout liste des distances dans le dictionnaire ?
#######  AUTRE PROBLÈME ==> index retourne la première occurence d'un élément et pas la liste des occurences
#######

def get_DistancePara(common_leaves, paralogs, matrix, dic):
    sum_para = 0
    i = -1
    
    for para in paralogs:
        if para in dic["Specie"][i+1:]:
            para1 = para +" "+dic["Identifier"][dic["Specie"].index(para, i+1, len(paralogs))]
            col = matrix[para1]
            for cl in common_leaves:
                leaf2 = cl+" "+dic["Identifier"][dic["Specie"].index(cl)]
                sum_para += col[leaf2]
            i = dic["Specie"].index(para, i+1, len(paralogs))
        else : i = -1
'''
