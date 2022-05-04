import holoviews as hv
import param
import panel as pn
from holoviews import opts, dim
import pandas as pd
import numpy as np
import seaborn as sns
from io import StringIO
from bokeh.models import HoverTool

from draggable_graph import DraggableGraph
from data_filter import DataFilter
from legends import nodes_colorbar
from utils import scale

class Network(param.Parameterized):
    
    node_data = param.DataFrame(precedence=-1)
    edge_data = param.DataFrame(precedence=-1) # list of [nodes, edges]
    graph_opts = param.Parameter({}, precedence=-1) # dict of {k: v} where k is a valid opts.k method and v is a dict of options
    selected_node = param.Tuple((None, None)) # tuple of (index_col, label_col) values for a selected node
    sel_nodes = param.DataFrame(precedence=-1)
    
    # data streams push data to DynamicMaps
    network_data = param.ClassSelector(default=hv.streams.Pipe(), class_=(hv.streams.Pipe,), precedence=-1) #ultimately this will be a list of node_data, edge_data, layout, bundle_edge_graphs
    click_stream = param.ClassSelector(default=hv.streams.Tap(), class_=(hv.streams.Tap,), precedence=-1)
    
    # graph layout algorithm
    layout = param.Selector(objects = ['kamada_kawai', 'circular', 'spring'], default='kamada_kawai')
    
    # edge bundling
    bundle_graph_edges = param.Selector(objects = ['Yes', 'No'], default='No')
    
    # network styling
    fontsize = param.Selector(default = '10pt', objects = ['{}pt'.format(i) for i in range(31)])
    node_cmap = param.Selector(objects = ['Category10', 
                                          'Category20', 
                                          'colorblind',
                                          'Holoviews', 
                                          'Set2', 
                                          'Set3', 
                                          'HTT_OMNI', 
                                          'YlGnBu', 
                                          'YlOrRed', 
                                          'Bokeh', 
                                          'isolum', 
                                          'Spectral', 
                                          'RdBu_r', 
                                          'RdGy_r'], 
                               default='HTT_OMNI')
    node_color = param.Selector(default = 'connectivity', objects = ['connectivity', 'data_source'])
    tooltips = param.ListSelector(default = [])
    label_color = param.Selector(default = 'black', objects = ['black', 'white'])
    cmap_centered = param.Boolean(default=False)
    node_clim = param.Tuple(default  = (None, None), precedence=-1)
    min_node_size = param.Number(default=25, bounds = (0, 150))
    max_node_size = param.Number(default=60, bounds = (0, 150))
    
    # parent DataFiter
    parent = param.ClassSelector(DataFilter, precedence=-1)
    
    loading = param.Boolean(default=False)
    click_loading = param.Boolean(default=False)
    
    def __init__(self, 
                 nodes, 
                 edges,
                 index_col = 'GeneID',
                 source_col = 'GENE_ID_A',
                 target_col = 'GENE_ID_B',
                 label_col = 'geneSymbol',
                 user_tooltips = [], # list of tuples (label, @column)
                 **params
                ):
        super(Network, self).__init__(**params)
        
        # TO DO: Ensure that any style parameters passed in graph.opts are equivalently passed as params during init
        
        self.nodes = nodes
        self.edges = edges
        
        self.index_col = index_col
        self.source_col = source_col
        self.target_col = target_col
        self.label_col = label_col
        
        # widget mappings
        self.mapping = dict([
            ('bundle_graph_edges', {'type': pn.widgets.RadioButtonGroup}
            ),
            ('fontsize', {'type': pn.widgets.DiscreteSlider, 'throttled': True}
            ),
            ('tooltips', {'type': pn.widgets.CheckBoxGroup}
            ),
            ('min_node_size', {'type': pn.widgets.FloatSlider, 'throttled': True, 'step': 5}
            ),
            ('max_node_size', {'type': pn.widgets.FloatSlider, 'throttled': True, 'step': 5}
            ),
        ]) 
        
        # export to csv
        self.export_show_nodes_button = pn.widgets.FileDownload(callback = self.export_show_nodes, filename = 'current_network_nodes.tab', label= 'Export current network nodes', button_type = 'primary')
        self.export_show_edges_button = pn.widgets.FileDownload(callback = self.export_show_edges, filename = 'current_network_edges.tab', label= 'Export current network edges', button_type = 'primary')
        self.export_sel_nodes_button = pn.widgets.FileDownload(callback = self.export_sel_nodes, filename = 'all_filtered_nodes.tab', label= 'Export all filtered nodes', button_type = 'default')
        self.export_sel_edges_button = pn.widgets.FileDownload(callback = self.export_sel_edges, filename = 'all_filtered_edges.tab', label= 'Export all filtered edges', button_type = 'default')

        self.network_pane = pn.pane.HoloViews(sizing_mode = 'stretch_both', linked_axes=False, min_height=0, min_width=0)
        self.cbar_pane = pn.pane.HoloViews(sizing_mode='stretch_width', linked_axes=False, min_height=0, min_width=0)
        
        # configure tooltips
        self.param.tooltips.objects = [i[0] for i in user_tooltips]+[self.parent.filter_aliases[f] for f in self.parent.filters]
        self.tooltip_map = {i:j for i, j in user_tooltips}
        self.tooltip_map.update({self.parent.filter_aliases[f]: f'@{f}' for f in self.parent.filters})
        
        self.graph = DraggableGraph(
            index_col = self.index_col,
            source_col = self.source_col,
            target_col = self.target_col,
            label_col = self.label_col,
        )
        
        ### configure cmap & node size ###
        if self.node_cmap == 'HTT_OMNI':
            node_cmap = sns.blend_palette(['white', '#4489ab'], as_cmap=True)
        elif self.node_cmap == 'Holoviews':
            node_cmap = hv.Cycle.default_cycles["default_colors"]
        else:
            node_cmap = self.node_cmap
        
        if not 'Nodes' in self.graph_opts:
            self.graph_opts['Nodes'] = {}
        self.graph_opts['Nodes'].update({'cmap': node_cmap})
        if not 'Graph' in self.graph_opts:
            self.graph_opts['Graph'] = {}
        self.graph_opts['Graph'].update({'cmap': node_cmap})
        
        self.param.node_color.objects = self.parent.color_opts

        max_PPI = 20
        size_dict = dict(zip(range(1, max_PPI+1), np.linspace(self.min_node_size, self.max_node_size, max_PPI)))
        size_dict.update(dict(zip(range(max_PPI+1, 101), [self.max_node_size]*len(list(range(max_PPI+1, 101))))))

        if not 'Nodes' in self.graph_opts:
            self.graph_opts['Nodes'] = {}
        self.graph_opts['Nodes'].update({'size': dim('PPI_SUM_TOTAL').categorize(size_dict)})
        if not 'Graph' in self.graph_opts:
            self.graph_opts['Graph'] = {}
        self.graph_opts['Graph'].update({'node_size': dim('PPI_SUM_TOTAL').categorize(size_dict)})

        # set up ordered param watching for cmap_centered, node
        self.param.watch(self.set_graph_opts_cmap_centered, 'cmap_centered', queued=True, precedence=2)
        
        self.update_sel_nodes()
        self.update_nodes_edges() # triggers self.update_data
        self.make_network_cbar()
        self.center_clim_bounds()
        
        # set tooltips to user provided tooltips
        self.tooltips = [i[0] for i in user_tooltips] # triggers self.update_tooltips
                 
    @param.depends('network_data.data', 'graph_opts', watch=True) # this method is the money maker!
    def view(self):
        self.loading=True
        network_graph = self.graph.view(self.network_data.data)
        
        if network_graph is not None:
            self.click_stream.source = self.graph.edge_graph
            network_graph = network_graph.opts(*[getattr(opts, k)(**self.graph_opts[k]) for k in self.graph_opts]).opts(active_tools=['point_draw'])

            if self.selected_node == (None, None): # only on initialization
                self.selected_node = tuple(self.node_data.iloc[0,:][[self.index_col, self.label_col]])
        
        self.network_pane.object = network_graph

        self.loading = False
        
    @param.depends('node_data', 'edge_data', 'layout', 'bundle_graph_edges', watch=True)
    def update_data(self):
        
        # for loading spinner control
        self.loading = True
        
        new_data = [
            self.node_data, # nodes
            self.edge_data, # edges
            self.layout,
            {'Yes': True, 'No': False}[self.bundle_graph_edges]            
        ]
        
        self.network_data.update(data=new_data) # triggers self.view
            
    @param.depends('click_stream.x', 'click_stream.y', watch=True)
    def get_clicked_node(self):
        
        self.click_loading = True
        
        x, y = (self.click_stream.x, self.click_stream.y)
        graph_data = pd.DataFrame(self.graph.stream.data)
        
        sel = graph_data.loc[((graph_data.x-x).abs()+(graph_data.y-y).abs()).idxmin(), :]
        idx = sel[self.index_col]
        label = sel[self.label_col]
        self.delta = (graph_data.x-x).abs()+(graph_data.y-y).abs().min()
        self.selected_node = (idx, label)
        
        self.click_loading = False
        
    @param.depends('parent.show_nodes', 'parent.show_edges', watch=True)
    def update_nodes_edges(self):
                
        new_nodes = self.parent.show_nodes.copy()
        new_edges = self.parent.show_edges.copy()
        
        new_edges['edge_width'] = scale(new_edges[self.parent.edge_score_col], 0.25, 5)
                
        self.param.set_param(node_data = new_nodes, edge_data = new_edges) # triggers self.update_data
        
    @param.depends('parent.sel_nodes', watch=True)
    def update_sel_nodes(self):
        self.param.set_param(sel_nodes = self.parent.sel_nodes.reset_index())
        
    @param.depends('fontsize', watch=True)
    def update_fontsize(self):
        # for loading spinner control
        self.loading = True
        
        if not 'Labels' in self.graph_opts:
            self.graph_opts['Labels'] = {}
        self.graph_opts['Labels'].update({'text_font_size': self.fontsize})
        self.param.set_param(graph_opts = self.graph_opts)
    
    @param.depends('node_cmap', watch=True)
    def update_node_cmap(self):
        # for loading spinner control
        self.loading = True
        
        if self.node_cmap == 'HTT_OMNI':
            node_cmap = sns.blend_palette(['white', '#4489ab'], as_cmap=True)
        elif self.node_cmap == 'Holoviews':
            node_cmap = hv.Cycle.default_cycles["default_colors"]
        else:
            node_cmap = self.node_cmap
        
        if not 'Nodes' in self.graph_opts:
            self.graph_opts['Nodes'] = {}
        self.graph_opts['Nodes'].update({'cmap': node_cmap})
        if not 'Graph' in self.graph_opts:
            self.graph_opts['Graph'] = {}
        self.graph_opts['Graph'].update({'cmap': node_cmap})
        
        self.param.set_param(graph_opts = self.graph_opts)
        
    @param.depends('node_color', watch=True)
    def update_node_color(self): 
        
        self.loading = True
        
        if not 'Nodes' in self.graph_opts:
            self.graph_opts['Nodes'] = {}
        self.graph_opts['Nodes'].update({'color': self.node_color})
        if not 'Graph' in self.graph_opts:
            self.graph_opts['Graph'] = {}
        self.graph_opts['Graph'].update({'node_color': self.node_color})
        
        self.param.set_param(graph_opts = self.graph_opts)
        
    @param.depends('node_clim', watch=True)
    def update_node_clim(self):
        
        if not 'Nodes' in self.graph_opts:
            self.graph_opts['Nodes'] = {}
        self.graph_opts['Nodes'].update({'clim': self.node_clim})
        if not 'Graph' in self.graph_opts:
            self.graph_opts['Graph'] = {}
        self.graph_opts['Graph'].update({'clim': self.node_clim})

    def set_graph_opts_cmap_centered(self, events):
        self.param.set_param(graph_opts = self.graph_opts)
        
    @param.depends('cmap_centered', 'node_color', 'parent.show_nodes', watch = True)
    def center_clim_bounds(self):
        
        if pd.api.types.is_numeric_dtype(self.parent.show_nodes[self.node_color]):
            min_, max_ = (self.parent.show_nodes[self.node_color].min(), self.parent.show_nodes[self.node_color].max())
            if self.cmap_centered == False:
                clim = (min_, max_)
            else:
                if (self.parent.show_nodes[self.node_color]>0).any() and (self.parent.show_nodes[self.node_color]<0).any():
                    bound = max(abs(min_), abs(max_))
                    clim = (-bound, bound)
                else:
                    clim = (min_, max_)
            
            self.node_clim = clim
        
        else:
            self.node_clim = (None, None)
         
    @param.depends('label_color', watch=True)
    def update_label_color(self):
        
        self.loading=True
        
        if not 'Labels' in self.graph_opts:
            self.graph_opts['Labels'] = {}
        self.graph_opts['Labels'].update({'text_color': self.label_color})
        
        self.param.set_param(graph_opts = self.graph_opts)
        
    @param.depends('parent.loading', watch=True)
    def update_loading(self):
        if self.parent.loading == True:
            self.loading = True
            
    @param.depends('loading', watch=True)
    def network_loading(self):
        self.network_pane.loading = self.loading
        
    @param.depends('node_color', 'node_cmap', 'node_clim', watch = True)
    def make_network_cbar(self):
        
        # show colorbar only if dtype of node_color column is numeric
        if pd.api.types.is_numeric_dtype(self.parent.show_nodes[self.node_color]):
            cbar = nodes_colorbar(cmap=self.graph_opts['Nodes']['cmap'], clim = self.node_clim, by=self.node_color)
            self.cbar_pane.object = cbar
            self.cbar_pane.visible = True
        else:
            self.cbar_pane.visible = False          
    
    @param.depends('tooltips', watch=True)
    def update_tooltips(self):
        
        hover = HoverTool(tooltips=[(i, self.tooltip_map[i]) for i in self.tooltips])
        
        if not 'Nodes' in self.graph_opts:
            self.graph_opts['Nodes'] = {}
        self.graph_opts['Nodes'].update({'tools': [hover]})
        
        self.param.set_param(graph_opts = self.graph_opts)
        
    @param.depends('parent.color_opts', watch = True)
    def update_color_opts(self):
        if not self.node_color in self.parent.color_opts:
            self.node_color = 'connectivity'
        self.param.node_color.objects = self.parent.color_opts
        
    @param.depends('min_node_size', 'max_node_size', watch = True)
    def update_node_size(self):
        max_PPI = 20
        size_dict = dict(zip(range(1, max_PPI+1), np.linspace(self.min_node_size, self.max_node_size, max_PPI)))
        size_dict.update(dict(zip(range(max_PPI+1, 101), [self.max_node_size]*len(list(range(max_PPI+1, 101))))))

        if not 'Nodes' in self.graph_opts:
            self.graph_opts['Nodes'] = {}
        self.graph_opts['Nodes'].update({'size': dim('PPI_SUM_TOTAL').categorize(size_dict)})
        if not 'Graph' in self.graph_opts:
            self.graph_opts['Graph'] = {}
        self.graph_opts['Graph'].update({'node_size': dim('PPI_SUM_TOTAL').categorize(size_dict)})

        self.param.set_param(graph_opts = self.graph_opts)

    def export_show_nodes(self):
        sio = StringIO()
        self.parent.show_nodes.to_csv(sio, sep='\t')
        sio.seek(0)
        
        return sio
    
    def export_show_edges(self):
        sio = StringIO()
        self.parent.show_edges.to_csv(sio, sep='\t')
        sio.seek(0)
        
        return sio
    
    def export_sel_nodes(self):
        sio = StringIO()
        self.parent.sel_nodes.to_csv(sio, sep='\t')
        sio.seek(0)
        
        return sio
    
    def export_sel_edges(self):
        sio = StringIO()
        self.parent.sel_edges.to_csv(sio, sep='\t')
        sio.seek(0)
        
        return sio