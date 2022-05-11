import numpy as np
import pandas as pd

def scale(arr, mn, mx, arr_min = None, arr_max = None):
    if len(arr)>1:
        arr = np.array(arr)

        if arr_min is None:
            arr_min = arr.min()
        else:
            arr = np.where(arr<arr_min, arr_min, arr)

        if arr_max is None:
            arr_max = arr.max()
        else:
            arr = np.where(arr>arr_max, arr_max, arr)

        return (mx-mn)*((arr-arr_min)/(arr_max-arr_min))+mn
    elif len(arr)==1:
        return np.array([mx])
    else:
        return arr

def save_hook(plot, element):
    plot.state.output_backend = 'svg'

def update_STRINGdb_edgefile(aliases_fn, links_fn):
    '''
    string_edgefile = update_STRINGdb_edgefile('9606.protein.aliases.v11.5.txt.gz', '9606.protein.links.v11.0.txt.gz')
    string_edgefile.to_csv('STRINGdb_edgefile.csv.gz', index=False)
    
    takes ~2 min to run -> STRINGdb_edgefile.csv.gz (~80 MB; unzipped -> ~500 MB)
    
    '''

    aliases = pd.read_csv(aliases_fn, sep='\t').dropna()

    temp = aliases[aliases['alias'].str.isnumeric()].groupby(['alias', '#string_protein_id']).size()
    strID_to_gID = temp.groupby('#string_protein_id').apply(lambda x: x.index.get_level_values(0).str.cat(sep=';'))
    
    interactions = pd.read_csv(links_fn, sep=' ', compression='gzip', index_col=[0,1]).squeeze().sort_index()
    interactions.index.names = ['stringId_A', 'stringId_B']
    interactions = pd.concat([interactions/1000, strID_to_gID.reindex(interactions.index, level='stringId_A'), strID_to_gID.reindex(interactions.index, level='stringId_B')], axis=1, keys = ['combined_score', 'GENE_ID_A', 'GENE_ID_B']).dropna()
    
    a = interactions['GENE_ID_A'].str.contains(';')
    b = interactions['GENE_ID_B'].str.contains(';')

    a_only = interactions[a&~b].groupby('stringId_A', group_keys=False).apply(lambda x: x.reset_index().set_index(x.index.names+['GENE_ID_B', 'combined_score'])['GENE_ID_A'].str.split(';', expand=True).stack())
    a_only.name = 'GENE_ID_A'

    b_only = interactions[b&~a].groupby('stringId_B', group_keys=False).apply(lambda x: x.reset_index().set_index(x.index.names+['GENE_ID_A', 'combined_score'])['GENE_ID_B'].str.split(';', expand=True).stack())
    b_only.name = 'GENE_ID_B'

    a_and_b = interactions[a&b].groupby('stringId_A', group_keys=False).apply(lambda x: x.reset_index().set_index(x.index.names+['GENE_ID_B', 'combined_score'])['GENE_ID_A'].str.split(';', expand=True).stack())
    a_and_b.name = 'GENE_ID_A'
    a_and_b = a_and_b.reset_index().set_index(interactions.index.names).groupby('stringId_B', group_keys=False).apply(lambda x: x.reset_index().set_index(x.index.names+['GENE_ID_A', 'combined_score'])['GENE_ID_B'].str.split(';', expand=True).stack())
    a_and_b.name = 'GENE_ID_B'
    
    edges = pd.concat([interactions[~(a|b)].reset_index(), a_only.reset_index(), b_only.reset_index(), a_and_b.reset_index()]).drop('level_4', axis=1)
    edges['EDGE_ID'] = np.where(edges['GENE_ID_A'].astype(float).astype(int)<edges['GENE_ID_B'].astype(float).astype(int), edges['GENE_ID_A'].str.cat(edges['GENE_ID_B'], sep=';'), edges['GENE_ID_B'].str.cat(edges['GENE_ID_A'], sep=';'))
    edges['GENE_ID_A'] = edges['GENE_ID_A'].astype(int)
    edges['GENE_ID_B'] = edges['GENE_ID_B'].astype(int)
    edges = edges[~edges['EDGE_ID'].duplicated()]

    # get rid of self edges
    edges = edges[edges['GENE_ID_A']!=edges['GENE_ID_B']]
    
    return edges