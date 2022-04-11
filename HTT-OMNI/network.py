import holoviews as hv
import param
import panel as pn
from holoviews import opts
import pandas as pd
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
    node_cmap = param.Selector(objects = hv.plotting.util.list_cmaps()+['HTT_OMNI'], default='HTT_OMNI')
    node_clim = param.Tuple(default = (None, None))
    tooltips = param.ListSelector(default = [])
    
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
        ]) 
        
        # export to csv
        self.export_show_nodes_button = pn.widgets.FileDownload(callback = self.export_show_nodes, filename = 'current_network_nodes.tab', label= 'Export current network nodes')
        self.export_show_edges_button = pn.widgets.FileDownload(callback = self.export_show_edges, filename = 'current_network_edges.tab', label= 'Export current network edges')
        self.export_sel_nodes_button = pn.widgets.FileDownload(callback = self.export_sel_nodes, filename = 'all_filtered_nodes.tab', label= 'Export all filtered nodes')
        self.export_sel_edges_button = pn.widgets.FileDownload(callback = self.export_sel_edges, filename = 'all_filtered_edges.tab', label= 'Export all filtered edges')

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
        
        if self.node_cmap == 'HTT_OMNI':
            node_cmap = sns.blend_palette(['white', '#4489ab'], as_cmap=True)
        else:
            node_cmap = self.node_cmap
        
        if not 'Nodes' in self.graph_opts:
            self.graph_opts['Nodes'] = {}
        self.graph_opts['Nodes'].update({'cmap': node_cmap})
        if not 'Graph' in self.graph_opts:
            self.graph_opts['Graph'] = {}
        self.graph_opts['Graph'].update({'cmap': node_cmap})
        
        self.update_sel_nodes()
        self.update_nodes_edges() # triggers self.update_data
        
        # set tooltips to user provided tooltips
        self.tooltips = [i[0] for i in user_tooltips] # triggers self.update_tooltips
                 
    @param.depends('network_data.data', 'graph_opts', watch=True) # this method is the money maker!
    def view(self):
        self.loading=True
        network_graph = self.graph.view(self.network_data.data).opts(*[getattr(opts, k)(**self.graph_opts[k]) for k in self.graph_opts]).opts(active_tools=['point_draw'])
        
        if network_graph is not None:
            self.click_stream.source = self.graph.edge_graph

            if self.selected_node == (None, None): # only on initialization
                self.selected_node = tuple(self.node_data.iloc[0,:][[self.index_col, self.label_col]])
        
        self.network_pane.object = network_graph
#         self.update_graph_opts() # Updates opts and sets self.loading to False
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
    
#     @param.depends('graph_opts', watch=True)
#     def update_graph_opts(self):
        
#         self.loading = True
        
#         if self.network_graph is not None:
#             self.network_pane.object = self.network_pane.object.opts(*[getattr(opts, k)(**self.graph_opts[k]) for k in self.graph_opts]).opts(active_tools=['point_draw'])
        
#         self.loading = False
        
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
        
        self.node_clim = (new_nodes['connectivity'].min(), new_nodes['connectivity'].max()) # triggers self.update_node_clim
        
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
        else:
            node_cmap = self.node_cmap
        
        if not 'Nodes' in self.graph_opts:
            self.graph_opts['Nodes'] = {}
        self.graph_opts['Nodes'].update({'cmap': node_cmap})
        if not 'Graph' in self.graph_opts:
            self.graph_opts['Graph'] = {}
        self.graph_opts['Graph'].update({'cmap': node_cmap})
        
        self.param.set_param(graph_opts = self.graph_opts)
        
    @param.depends('node_clim', watch=True)
    def update_node_clim(self): # Note: graph_opts must be set after call (e.g. after any change to self.node_clim)
        
        if not 'Nodes' in self.graph_opts:
            self.graph_opts['Nodes'] = {}
        self.graph_opts['Nodes'].update({'clim': self.node_clim})
        if not 'Graph' in self.graph_opts:
            self.graph_opts['Graph'] = {}
        self.graph_opts['Graph'].update({'clim': self.node_clim})
        
    @param.depends('parent.loading', watch=True)
    def update_loading(self):
        if self.parent.loading == True:
            self.loading = True
            
    @param.depends('loading', watch=True)
    def network_loading(self):
        self.network_pane.loading = self.loading
        
    @param.depends('node_clim', 'node_cmap', watch = True)
    def make_network_cbar(self):
        cbar = nodes_colorbar(cmap=self.graph_opts['Nodes']['cmap'], clim = self.node_clim)
        self.cbar_pane.object = cbar
    
    @param.depends('tooltips', watch=True)
    def update_tooltips(self):
        
        hover = HoverTool(tooltips=[(i, self.tooltip_map[i]) for i in self.tooltips])
        
        if not 'Nodes' in self.graph_opts:
            self.graph_opts['Nodes'] = {}
        self.graph_opts['Nodes'].update({'tools': [hover]})
        
        self.param.set_param(graph_opts = self.graph_opts)
        
    ##### need to configure these to only include relevant columns ######
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