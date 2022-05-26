import pandas as pd
import numpy as np
import holoviews as hv
from holoviews import opts, dim
import panel as pn
import os
import dask.dataframe as dd

def setup():
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
    pn.state.cache = {}

    def unpack_STRINGdb_edgefile():
        stringdb_edgefile = pd.read_csv(r'./assets/data/STRINGdb_edgefile.csv.gz')
        stringdb_edgefile.to_csv(r'./assets/data/STRINGdb_edgefile.csv', index=False)

    # since STRINGdb_edgefile is too large to track using normal git, 
    # we'll just git track the gzipped version and unpack it locally when needed
    if not os.path.exists(r'./assets/data/STRINGdb_edgefile.csv'):
        unpack_STRINGdb_edgefile()

    def save_hook(plot, element):
        plot.state.output_backend = 'svg'

    geneID_col = 'interactor_Human_Ortholog_EntrezGeneID'
    geneSymbol_col = 'interactor_Human_Ortholog_EntrezGeneSymbol'
    
   ################################ READ IN NODES ################################

    nodes = pd.read_csv(r'./assets/data/nodes.csv', low_memory=False)
    nodes.columns = nodes.columns.str.replace(' ', '_')

    # remove rows with only a negative interaction result
    nodes = nodes[nodes['interaction_result'].str.contains('Y').fillna(True)]

    # replace Cell Culture with Cell culture
    nodes['model'] = nodes['model'].str.replace('Cell Culture', 'Cell culture')

    # replace "brain" tissue annotation with "brain, whole"
    nodes['tissue'] = nodes['tissue'].where(nodes['tissue']!='brain', 'brain (whole)')

    # combine model and model_organism annotations (for Cell culture and In vitro)
    nodes['model_species'] = nodes['model'].where(~nodes['model'].isin(['Cell culture', 'In vitro']), nodes['model'].str.cat(nodes['model_organism'].fillna('Not reported'), sep=' (')+')')

    # create additional annotation for AP-MS studies
    nodes['detection_method_annot'] = nodes['base_method'].where(lambda x: nodes['base_method']!='affinity chromatography technology', nodes['base_method'].str.cat(nodes['detection_method'], sep=' (')+')') 

    # drop any rows without a human ortholog
    nodes = nodes[~nodes[geneID_col].isnull()]

    # some dtype mapping
    nodes[geneID_col] = nodes[geneID_col].astype(int)
    nodes['year'] = nodes['year'].astype(int)

    # consolidate multiple GeneSymbols for a single GeneID
    newest_geneSymbols = nodes.groupby(geneID_col).apply(lambda x: x.loc[x['year'].idxmax(), geneSymbol_col])
    nodes[geneSymbol_col] = nodes[geneID_col].map(newest_geneSymbols)

    # fill in common_name column with "WT" for PubMed:22556411
    nodes['common_name'] = nodes['common_name'].where(nodes['source_identifier']!='PubMed:22556411', 'WT')

    # add study_id column
    nodes['study_id'] = nodes['authors'].str.split(',', n=1, expand=True)[0].str.split(' ', expand=True, n=1)[1].str.cat(nodes['year'].astype(str), sep=' ').str.cat(nodes['journal'], sep=' ').fillna(nodes['source_identifier'])
    
    ################################ READ IN OMICS DATA ################################

    omics_data = pd.read_csv(r'./assets/data/20220319_omics_data.csv', header=[0, 1, 2, 3], index_col=[0, 1])
    omics_data.index.names = ['geneSymbol', 'geneID']
    omics_data = omics_data.reset_index().set_index(['geneID', 'geneSymbol'])
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

    background_geneIDs = omics_data.index.get_level_values('geneID').astype(str).values.tolist()

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
    
    # fill in missing filter values with "Not reported"
    nodes['data_source'] = 'HINT'
    nodes.loc[:, filters] = nodes.loc[:, filters].where(lambda x: x.notnull(), 'Not reported')

    ################################# NETWORK ##############################
    node_color = 'connectivity'
    node_cmap = 'HTT_OMNI'
    tooltips = [
        ('Gene ID', f'@{geneID_col}'),
        ('Gene Symbol', f'@{geneSymbol_col}'), 
        ('# PPI observations (all)', '@PPI_SUM_TOTAL'),
    ]

    graph_opts = {
        'Nodes': dict(
            color=node_color,
            cmap = node_cmap,
        ),
        'Graph': dict(
            edge_color = 'grey',
            edge_line_width = 'edge_width',
            node_color = node_color,
            cmap = node_cmap,
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
        'GO biological process': 'GO:0008150',
        'GO molecular function': 'GO:0003674',
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



    def tissue_age_leg(all_tissues, all_ages, x, graph_opts = []):
        '''
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

        '''

        tissue_age = all_tissues.tolist()+all_ages.tolist()
        
        dummy = pd.DataFrame(np.zeros(shape=(len(tissue_age), 2)), index=tissue_age, columns=[x, x])
        dummy_ds = hv.Dataset(dummy.reset_index().melt('index'))
        dummy_ds.columns = ['tissue_age', 'variable', 'value']

        return dummy_ds.to(hv.Curve, 'variable', ['value'], 'index').overlay('index').opts(*graph_opts)

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

    # cache required data
    pn.state.cache['nodes'] = nodes
    pn.state.cache['edges'] = dd.read_csv(r'./assets/data/STRINGdb_edgefile.csv', usecols=[2, 3, 4])
    pn.state.cache['filters'] = filters
    pn.state.cache['index_col'] = geneID_col
    pn.state.cache['gene_symbol_col'] = geneSymbol_col
    pn.state.cache['filter_aliases'] = filter_aliases
    pn.state.cache['groupby_PPI_cols'] = [geneID_col, 'source_identifier']

    pn.state.cache['graph_opts'] = graph_opts
    pn.state.cache['source_col'] = 'GENE_ID_A'
    pn.state.cache['target_col'] = 'GENE_ID_B'
    pn.state.cache['label_col'] = geneSymbol_col
    pn.state.cache['fontsize'] = graph_opts['Labels']['text_font_size']
    pn.state.cache['node_cmap'] = node_cmap
    pn.state.cache['user_tooltips'] = tooltips

    pn.state.cache['annot_description_mapping'] = annot_desc

    pn.state.cache['omics_data'] = omics_data
    pn.state.cache['dummy_leg'] = dummy_leg
    pn.state.cache['plot_opts'] = plot_opts
    pn.state.cache['background_geneIDs'] = background_geneIDs

    if os.path.exists(r'./assets/data/init_GO_results.csv'):
        pn.state.cache['GO_init_results'] = pd.read_csv(r'./assets/data/init_GO_results.csv')


setup()