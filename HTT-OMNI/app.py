import panel as pn
import param

from data_filter import DataFilter
from network import Network
from enrichment import Enrichment
from omics_data_viewer import OmicsDataViewer

class App(param.Parameterized):
    
    data_filter = param.ClassSelector(DataFilter)
    network = param.ClassSelector(Network)
    enrichment = param.ClassSelector(Enrichment)
    omics_viewer = param.ClassSelector(OmicsDataViewer)
    
    def __init__(self, **params):
        super(App, self).__init__(**params)

        COLS = 12
        ROW_HEIGHT = 75
        ACCENT_BASE_COLOR = '#4489ab'
        sidebar_width = 500

        react = pn.template.ReactTemplate(
            title='Explore HTT Interactions', 
            sidebar_width=sidebar_width, 
            row_height = ROW_HEIGHT, 
            prevent_collision = True,
            header_background = ACCENT_BASE_COLOR, 
            logo = r'.\assets\images\CHDI_logo_white.png',
            favicon = r'.\assets\images\CHDI_logo.png',
            cols={'lg': COLS, 'md': COLS, 'sm': COLS, 'xs': COLS, 'xxs': COLS},
        )
        
        widgets = {}

        widgets.update(self.enrichment.mapping.copy())
        widgets.update(self.data_filter.mapping.copy())
        widgets.update(self.network.mapping.copy())
        widgets.update(self.omics_viewer.mapping.copy())

        param_opts = {'widgets': widgets, 'show_name': False,}

        node_properties = [
            pn.Row(
                '**Show unconnected nodes**', 
                pn.Param(self.data_filter, parameters = ['vis_unconnected'], **param_opts), 
                height=15
            ), 
            pn.Param(self.data_filter, parameters = ['max_nodes', 'node_display_priority'], **param_opts),
        ]

        edge_properties = [
            pn.Param(self.data_filter, parameters = ['STRINGdb_score'], **param_opts), 
            pn.Param(self.network, parameters = ['layout'], **param_opts),
            pn.Row(
                pn.pane.Markdown('Bundle edges?', align='center'), 
                pn.Param(self.network, parameters = ['bundle_graph_edges'], **param_opts),
                height = 15
            ),
            pn.Column()
        ]
        
        data_filters = [
            pn.Card(
                pn.Param(self.data_filter, parameters = [f+'_AND_OR_NOT', f], show_labels=False, **param_opts), 
                collapsed=True, 
                title=self.data_filter.mapping[f]['name'], 
            ) 
            for f in self.data_filter.filters]
        
        node_filter_wids = pn.Column(
            pn.Param(self.data_filter, parameters = ['PPI_sum_cutoff'], **param_opts),
            pn.Card(
                pn.Param(self.data_filter, parameters=['node_query'], show_labels=False, **param_opts), 
                collapsed=True, 
                title='Search for Gene IDs',
            ),
            *data_filters
        )
        
        aesthetic_properties = [
            pn.Param(self.network, parameters = ['fontsize', 'node_color', 'node_cmap', 'label_color', 'cmap_centered'], **param_opts),
            pn.Card(pn.Param(self.network, parameters=['tooltips'], **param_opts), title='Hover tooltip info')
        ]
        
        user_upload = [
            pn.Column(
                pn.Row(
                    'Select a file to upload', 
                    pn.Param(self.data_filter, parameters = ['user_upload_file'], **param_opts)
                ),
                'Required column headers: gene_id, string_id, study_id',
                'Optional column headers: '+', '.join([i for i in self.data_filter.filters if not i=='study_id']),
                
            )
        ]

        AS_toggles = pn.Card(pn.Column(pn.Param(self.omics_viewer, parameters = ['tissues', 'ages'], **param_opts)), title='AS omics data toggles')

        sidebar = pn.Tabs(
            pn.Column(
                pn.Card(*node_properties, title = 'Node Properties'), 
                pn.Card(*edge_properties, title = 'Edge Properties'),
                pn.Card(*aesthetic_properties, title = 'Aesthetic Properties'),
                name = 'NETWORK PROPERTIES'
            ), 
            pn.Column(
                pn.Param(self.data_filter, parameters = ['reset_filters'], **param_opts),
                node_filter_wids, 
                name = 'DATA FILTERS'
            ),
            pn.Column(
                AS_toggles,
                name = 'PLOT FILTERS'                
            ),
            pn.Column(*user_upload, name = 'UPLOAD DATA'),
            max_width = sidebar_width-50
        )

        react.sidebar.append(sidebar)

        react.main[:12,:6] = pn.Tabs(
            pn.Column(
                pn.Param(self.data_filter, parameters = ['network_plot_title'], show_labels=False, widgets = widgets, show_name = False),
                self.network.cbar_pane,
                self.network.network_pane,
                pn.Row(
                    self.network.export_show_nodes_button, 
                    self.network.export_show_edges_button, 
                    self.network.export_sel_nodes_button, 
                    self.network.export_sel_edges_button,
                ),
                name = 'Network visualization'
            ),
            pn.Column(
                pn.Row(
                    pn.Param(self.enrichment, parameters = ['GO_annot'], show_labels = False, **param_opts), 
                    pn.Param(self.enrichment, parameters = ['run_GO_analysis'], show_labels = False, **param_opts)
                ),
                pn.Param(self.enrichment, parameters = ['GO_plot_title'], show_labels=False, **param_opts),
                self.enrichment.plot_pane,
                pn.Param(self.enrichment, parameters = ['GO_show', 'GO_max_FDR', 'GO_min_enrichment'], default_layout=pn.Row, **param_opts),
                name = 'GO/pathway enrichment analysis',
            )
        )

        react.main[:6, 6:] = pn.Card(
            pn.Tabs(
                pn.Column(self.omics_viewer.selected_node_model_obs_pane, name = 'PPI per model'),
                pn.Column(self.omics_viewer.selected_node_AS_RNA_pane, name = 'RNA levels'),
                pn.Column(self.omics_viewer.selected_node_AS_protein_pane, name = 'Protein levels'),
                pn.Column(self.omics_viewer.selected_node_snRNA_pane, name = 'snRNA levels'),
                pn.Column(self.omics_viewer.selected_node_scRNA_pane, name = 'scRNA levels'),
                dynamic = True,
            ), 
            title = 'Omics data (selected node)',
        )
        
        react.main[6:12, 6:] = pn.Card(
            pn.Tabs(
                pn.Column(self.omics_viewer.network_model_obs_pane, name = 'PPI per model'),
                pn.Column(self.omics_viewer.network_AS_RNA_pane, name = 'RNA levels'),
                pn.Column(self.omics_viewer.network_AS_protein_pane, name = 'Protein levels'),
                pn.Column(self.omics_viewer.network_snRNA_pane, name = 'snRNA levels'),
                pn.Column(self.omics_viewer.network_scRNA_pane, name = 'scRNA levels'),
                dynamic = True
            ), 
            title = 'Omics data (filtered nodes)',
        )
        react.main[12:18, :] = pn.Card(pn.Param(self.data_filter, parameters = ['display_nodes'], **param_opts), title = 'Network nodes table')

        self.template = react        
        
    def view(self):
        self.template.servable()