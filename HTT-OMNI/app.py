import holoviews as hv
import panel as pn
from holoviews import opts, dim
import pandas as pd
import numpy as np
from bokeh.models import HoverTool

from network import Network
from enrichment import Enrichment
from omics_data_viewer import OmicsDataViewer
from legends import nodes_colorbar



nodes = pd.read_csv('20220310_nodes.csv.gz', low_memory=False)
edges = pd.read_csv('20220310_edges.csv.gz')

geneID_col = 'interactor_Human Ortholog_EntrezGeneID'
geneSymbol_col = 'interactor_Human Ortholog_EntrezGeneSymbol'

test_nodes = nodes.groupby([geneID_col, geneSymbol_col]).size().reset_index()
test_nodes.columns = ['geneID', 'geneSymbol', 'PPI_SUM']
test_nodes = test_nodes.sort_values('PPI_SUM', ascending=False).iloc[:50, :]
test_edges = edges[edges['GENE_ID_A'].isin(test_nodes['geneID'])&edges['GENE_ID_B'].isin(test_nodes['geneID'])]

smallest, largest = (25, 60)
max_PPI = 20
size_dict = dict(zip(range(1, max_PPI+1), np.linspace(smallest, largest, max_PPI)))
size_dict.update(dict(zip(range(max_PPI+1, 101), [largest]*len(list(range(max_PPI+1, 101))))))

omics_data = pd.read_csv(r'Z:\2 Programming\10 CHDI\PPI visualization\24mo milestone\20220319_omics_data.csv', header=[0, 1, 2, 3], index_col=[0, 1])
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

################################# NETWORK ##############################
node_color = 'connectivity'
node_size = dim('PPI_SUM').categorize(size_dict)
node_cmap = 'Blues'
tools = [HoverTool(tooltips = [('Gene ID', '@geneID'), ('Gene Symbol', '@geneSymbol'), ('# PPI observations (all)', '@PPI_SUM')])]

graph_opts = {
    'Nodes': dict(
        color=node_color,
        cmap = node_cmap,
        size = node_size,
#         tools = tools+['tap', 'point_draw']
#         clim = clim,
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
#         tools = tools+['tap', 'point_draw']
#         active_tools = ['point_draw']
    ),
    'Labels': dict(text_font_size = '10pt')
}

app = Network(
    nodes = test_nodes, 
    edges = test_edges,
    graph_opts = graph_opts, 
    index_col = 'geneID', 
    source_col = 'GENE_ID_A',
    target_col = 'GENE_ID_B',
    label_col = 'geneSymbol',
    fontsize = graph_opts['Labels']['text_font_size'],
    node_cmap = graph_opts['Nodes']['cmap']
)
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

ignore_GO_terms = ['biological_process', 'cellular process']

enrichment = Enrichment(app, annot_desc)  

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
    ],
}
omics_viewer = OmicsDataViewer(omics_data, parent = app, plot_opts=plot_opts)

cbar = hv.DynamicMap(nodes_colorbar, streams = dict(zip(['cmap', 'clim'],[app.param.node_cmap, app.param.node_clim])))
cbar_pane = pn.pane.HoloViews(cbar, sizing_mode='stretch_width', linked_axes=False, min_height=0, min_width=0)

widgets = dict([
    ('max_nodes', {'type': pn.widgets.DiscreteSlider, 'throttled': True}
    ),
    ('bundle_graph_edges', {'type': pn.widgets.RadioButtonGroup}
    ),
    ('fontsize', {'type': pn.widgets.DiscreteSlider, 'throttled': True}
    ),
    ('run_GO_analysis', {'button_type': 'default'}
    ),
    ('GO_annot', {'type': pn.widgets.Select,}
    ),
    ('GO_show', {'type': pn.widgets.IntSlider, 
              'throttled': True, 
              'name': 'Displayed # of enriched terms'}
    ),
    ('GO_max_FDR', {'type': pn.widgets.FloatInput, 
                 'name': 'FDR cutoff', 
                 'step': 0.01}
    ),
    ('GO_min_enrichment', {'type': pn.widgets.FloatInput,
                        'name': 'Minimum fold enrichment', 
                        'step': 0.5}
    ),
    ('GO_plot_title', {'type': pn.widgets.StaticText}
    ),
    ('tissues', {'type': pn.widgets.MultiChoice, 'name': 'Mouse AS tissues'}
    ), 
    ('ages', {'type': pn.widgets.MultiChoice, 'name': 'Mouse AS ages (months)'}
    ),
    ('age_label', {'type': pn.widgets.StaticText}
    ),
    ('tissue_label', {'type': pn.widgets.StaticText}
    )
])
param_opts = {'widgets': widgets, 'show_name': False}

css = """
.bk.card button.bk.card-header .bk.card-header-row .bk .bk.bk-clearfix {
  font-size: 16px;
  text-align: left;
  padding-left: 5px;
}
"""
pn.extension(raw_css=[css])
pn.extension(sizing_mode = 'stretch_width')

COLS = 12
ROW_HEIGHT = 75
ACCENT_BASE_COLOR = '#4489ab'

react = pn.template.ReactTemplate(
    title='Explore HTT Interactions', 
    sidebar_width=400, 
    row_height = ROW_HEIGHT, 
    prevent_collision = True,
    header_background = ACCENT_BASE_COLOR, 
    logo = 'CHDI_logo_white.png',
    favicon = 'CHDI_logo.png',
    cols={'lg': COLS, 'md': COLS, 'sm': COLS, 'xs': COLS, 'xxs': COLS},
)

react.sidebar.append(pn.Param(app, parameters = ['max_nodes', 'layout', 'bundle_graph_edges', 'fontsize', 'node_cmap'], **param_opts))
react.sidebar.append(pn.Card(pn.Column(pn.Param(omics_viewer, parameters = ['tissues', 'ages'], **param_opts)), title='AS omics data toggles'))

react.main[:10,:6] = pn.Tabs(
    pn.Column(
                cbar_pane,
        app.network_pane,
        name = 'Network visualization'
    ),
    pn.Column(
        pn.Row(
            pn.Param(enrichment, parameters = ['GO_annot'], show_labels = False, **param_opts), 
            pn.Param(enrichment, parameters = ['run_GO_analysis'], show_labels = False, **param_opts)
        ),
        pn.Param(enrichment, parameters = ['GO_plot_title'], show_labels=False, **param_opts),
        enrichment.plot_pane,
        pn.Param(enrichment, parameters = ['GO_show', 'GO_max_FDR', 'GO_min_enrichment'],default_layout=pn.Row, **param_opts),
        name = 'GO/pathway enrichment analysis',
    )
)

react.main[:6, 6:] = pn.Card(
    pn.Tabs(
        pn.Column(omics_viewer.selected_node_AS_RNA_pane, name = 'RNA'),
        pn.Column(omics_viewer.selected_node_AS_protein_pane, name = 'Protein'),
#         **tab_params
    ), 
    title = 'Omics data (selected node)',
)
react.main[6:12, 6:] = pn.Card(
    pn.Tabs(
        pn.Column(omics_viewer.network_AS_RNA_pane, name = 'RNA'),
        pn.Column(omics_viewer.network_AS_protein_pane, name = 'Protein'),
#         **tab_params
    ), 
    title = 'Omics data (filtered nodes)',
)
react.show()