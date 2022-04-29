import pandas as pd
import numpy as np
import holoviews as hv
from holoviews import opts, dim
import panel as pn
import os

from data_filter import DataFilter
from network import Network
from enrichment import Enrichment
from omics_data_viewer import OmicsDataViewer
from app import App
from utils import save_hook
from legends import tissue_age_leg

hv.extension('bokeh')

css = """
.bk.card button.bk.card-header .bk.card-header-row .bk .bk.bk-clearfix {
    font-size: 16px;
    text-align: left;
    padding-left: 5px;
}

.bk-root .bk-btn-primary {
    color: #fff; 
    background-color: #4489ab;
}

.bk-root .bk-btn-primary:hover {
      background-color: #366c87;
}
input[type=file] {
    width: 100%;
    height: 100px;
    border-width: 3px !important;
    border-style: dashed !important;
    border-color: #9E9E9E !important;
    background: #EEEEEE;
    border-radius: 5px;
    text-align: center;
    margin: auto;
}
"""

pn.extension(
    sizing_mode='stretch_width',
    loading_spinner='dots', 
    loading_color='#4489ab',
    notifications = True,
)

pn.config.raw_css.append(css)
pn.param.ParamMethod.loading_indicator = True

# since STRINGdb_edgefile is too large to track using normal git, 
# we'll just git track the gzipped version and unpack it locally when needed
if not os.path.exists(r'.\assets\data\STRINGdb_edgefile.csv'):
    stringdb_edgefile = pd.read_csv(r'.\assets\data\STRINGdb_edgefile.csv.gz')
    stringdb_edgefile.to_csv(r'.\assets\data\STRINGdb_edgefile.csv', index=False)

nodes = pd.read_csv(r'.\assets\data\20220319_published_nodes.csv.gz', low_memory=False)
nodes.columns = nodes.columns.str.replace(' ', '_')

edges = pd.read_csv(r'.\assets\data\20220319_published_edges.csv.gz')
edges.columns = edges.columns.str.replace(' ', '_')

geneID_col = 'interactor_Human_Ortholog_EntrezGeneID'
geneSymbol_col = 'interactor_Human_Ortholog_EntrezGeneSymbol'

test_nodes = nodes.groupby([geneID_col, geneSymbol_col]).size().reset_index()
test_nodes.columns = ['geneID', 'geneSymbol', 'PPI_SUM']
test_nodes = test_nodes.sort_values('PPI_SUM', ascending=False).iloc[:50, :]
test_edges = edges[edges['GENE_ID_A'].isin(test_nodes['geneID'])&edges['GENE_ID_B'].isin(test_nodes['geneID'])]

smallest, largest = (25, 60)
max_PPI = 20
size_dict = dict(zip(range(1, max_PPI+1), np.linspace(smallest, largest, max_PPI)))
size_dict.update(dict(zip(range(max_PPI+1, 101), [largest]*len(list(range(max_PPI+1, 101))))))

omics_data = pd.read_csv(r'.\assets\data\20220319_omics_data.csv', header=[0, 1, 2, 3], index_col=[0, 1])
omics_data.index.names = ['geneSymbol', 'geneID']
omics_data = omics_data.reset_index().set_index(['geneID', 'geneSymbol'])
omics_data = omics_data[omics_data.index.get_level_values('geneID').isin(nodes[geneID_col])].copy()
temp = omics_data.T.reset_index()
temp['Q-length'] = temp['Q-length'].str.strip('Q').astype(np.int64)
temp['age'] = temp['age'].astype(np.int64)
temp['Tissue/Age'] = temp['tissue']+' ('+temp['age'].astype(str)+'mo)'
omics_data = temp.set_index(omics_data.columns.names+['Tissue/Age']).T

all_tissues = np.unique(omics_data[['PROTEIN', 'RNA']].columns.get_level_values('tissue'))
all_ages = np.unique(omics_data[['PROTEIN', 'RNA']].columns.get_level_values('age'))

tissue_colors = dict(zip(all_tissues, hv.Cycle.default_cycles["default_colors"]))
tissue_colors.update(dict(zip(all_ages, ['#000000']*len(all_ages))))

age_dashes = dict(zip(all_ages, ['dotted', 'dashed', 'solid']))
age_dashes.update(dict(zip(all_tissues, ['solid']*len(all_tissues))))

############################### DATA FILTER ############################
filters = [
    'model_species', 
    'common_name',
    'cell_culture_comment',
    'tissue',
    'htt_length',
    'detection_method_annot',
    'study_id',
    'data_source'
]

filter_aliases = dict(zip(filters, ['Model (species)', 'Mouse model ID', 'Cell culture subtype', 'Tissue', 'HTT length', 'Method', 'Study (first author, year, journal)', 'Data source']))

################################# NETWORK ##############################
node_color = 'connectivity'
node_size = dim('PPI_SUM_TOTAL').categorize(size_dict)
node_cmap = 'HTT_OMNI'
tooltips = [
    ('Gene ID', f'@{geneID_col}'),
    ('Gene Symbol', f'@{geneSymbol_col}'), 
    ('# PPI observations (all)', '@PPI_SUM_TOTAL'),
    ('Data soource', '@data_source')
]

graph_opts = {
    'Nodes': dict(
        color=node_color,
        cmap = node_cmap,
        size = node_size,
    ),
    'Graph': dict(
        edge_color = 'grey',
        edge_line_width = 'edge_width',
        node_color = node_color,
        cmap = node_cmap,
        node_size = node_size,
#         clim = clim,
        tools = ['tap']
    ),
    'Overlay': dict(
        xlim=(-1.3, 1.3),
        ylim=(-1.3, 1.3),
        responsive=True,
        xaxis=None, 
        yaxis=None,
        min_height=0,
        min_width=0,
        hooks = [save_hook]
    ),
    'Labels': dict(text_font_size = '10pt')
}

################################# ENRICHMENT ##############################
annot_desc = {
    'GO biological process': 'GO:0003674',
    'GO molecular function': 'GO:0008150',
    'GO cellular component': 'GO:0005575',
    'GO SLIM molecular function': 'ANNOT_TYPE_ID_PANTHER_GO_SLIM_MF',
    'GO SLIM biological process': 'ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP',
    'GO SLIM cellular component': 'ANNOT_TYPE_ID_PANTHER_GO_SLIM_CC',
    'Panther pathways': 'ANNOT_TYPE_ID_PANTHER_PATHWAY',
    'Reactome pathways': 'ANNOT_TYPE_ID_REACTOME_PATHWAY'
}

ignore_GO_terms = ['biological_process', 'cellular process', 'cellular_component']

################################# OMICS PLOTS ##############################
plot_opts = {
    'PROTEIN': [
        opts.Curve(
            line_color = dim('tissue').categorize(tissue_colors),
            line_dash=dim('age').categorize(age_dashes),
            show_legend=False,
            ylabel = 'log2FC abundance',
            framewise=True,
            responsive=True,
            xlim = (45, 180)
        ),
        opts.Scatter(
            line_color = dim('tissue').categorize(tissue_colors),
            show_legend=False,
            ylabel = 'log2FC abundance',
            size=7,
            color='white',
            framewise=True,
            responsive=True,
            xlim = (45, 180)
        ),
        opts.Overlay(
            hooks = [save_hook]
        )
    ],
    'RNA': [
        opts.Curve(
            line_color = dim('tissue').categorize(tissue_colors),
            line_dash=dim('age').categorize(age_dashes),
            show_legend=False,
            ylabel = 'log2FC abundance',
            framewise=True,
            responsive=True,
            xlim = (75, 180)
        ),
        opts.Scatter(
            line_color = dim('tissue').categorize(tissue_colors),
            show_legend=False,
            ylabel = 'log2FC abundance',
            size=7,
            color='white',
            framewise=True,
            responsive=True,
            xlim = (75, 180)
        ),
        opts.Overlay(
            hooks = [save_hook]
        )
    ],
    'SCRNA':[
        opts.Bars(
            fill_color = '#4489ab',
            alpha='alpha',
            responsive=True, 
            show_frame=False,
            xrotation=90,
            ylim = (0, 1.05),
            hooks = [save_hook],            
        )
    ],
    'SNRNA':[
        opts.Bars(
            color = (dim('log2FC abundance (Q175/WT)')<0).categorize({True:'#4489ab', False:'#ab4444'}),
            line_color = (dim('log2FC abundance (Q175/WT)')<0).categorize({True:'#4489ab', False:'#ab4444'}),
            fill_alpha = 'alpha',
            responsive = True,
            framewise=True,
            xrotation=90,
            show_frame=False,
            hooks = [save_hook]
        ),
    ],
    'model_obs': [
        opts.Bars(
            fill_color='#4489ab', 
            responsive=True, 
            show_frame=False,
            framewise=True,
            xlabel = '',
            hooks = [save_hook]
        ),
    ]
}

dummy_leg_opts = [
    opts.Curve(
        line_color = dim('index').categorize(tissue_colors),
        line_dash=dim('index').categorize(age_dashes),
        framewise=True,    
    ),
    opts.NdOverlay(
        legend_position='right', 
        show_legend=True,
        legend_opts = {'title':'Tissue/Age (mo)'},
        framewise=True,
        show_frame = False,
        hooks = [save_hook]
    )
]

dummy_leg = tissue_age_leg(all_tissues, all_ages, 150, graph_opts = dummy_leg_opts)

def user_instance():

    data_filter = DataFilter(
        nodes = nodes, 
        edges = edges, 
        filters = filters, 
        index_col=geneID_col,
        gene_symbol_col = geneSymbol_col,
        filter_aliases = filter_aliases,
        groupby_PPI_cols = [geneID_col, 'source_identifier']
    )

    network = Network(
        parent = data_filter,
        nodes = nodes, 
        edges = edges,
        graph_opts = graph_opts, 
        index_col = geneID_col, 
        source_col = 'GENE_ID_A',
        target_col = 'GENE_ID_B',
        label_col = geneSymbol_col,
        fontsize = graph_opts['Labels']['text_font_size'],
        node_cmap = graph_opts['Nodes']['cmap'],
        user_tooltips = tooltips,
    )

    enrichment = Enrichment(network, annot_desc, index_col=geneID_col)  

    omics_viewer = OmicsDataViewer(omics_data, dummy_leg = dummy_leg, parent = network, plot_opts=plot_opts)

    app = App(
        enrichment = enrichment, 
        data_filter=data_filter, 
        network=network, 
        omics_viewer=omics_viewer,
    )
    
    return app.view()

user_instance()