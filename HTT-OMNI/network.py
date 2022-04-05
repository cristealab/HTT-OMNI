import holoviews as hv
import param
import panel as pn
from holoviews import opts
import pandas as pd

from draggable_graph import DraggableGraph
from utils import scale

class Network(param.Parameterized):
    
    node_data = param.DataFrame()
    edge_data = param.DataFrame() # list of [nodes, edges]
    graph_opts = param.Parameter({}) # dict of {k: v} where k is a valid opts.k method and v is a dict of options
    sel_data = param.DataFrame()
    selected_node = param.Tuple((None, None)) # tuple of (index_col, label_col) values for a selected node
    
    # data streams push data to DynamicMaps
    network_data = param.ClassSelector(default=hv.streams.Pipe(), class_=(hv.streams.Pipe,), precedence=-1) #ultimately this will be a list of node_data, edge_data, graph_opts_data, layout, bundle_edge_graphs
    click_stream = param.ClassSelector(default=hv.streams.Tap(), class_=(hv.streams.Tap,), precedence=-1)
    
    # slider for # of nodes to show
    max_nodes = param.Selector(objects = [10, 20, 50, 100, 200, 300, 400, 500], default=50)
    
    # graph layout algorithm
    layout = param.Selector(objects = ['kamada_kawai', 'circular', 'spring'], default='kamada_kawai')
    
    # edge bundling
    bundle_graph_edges = param.Selector(objects = ['Yes', 'No'], default='No')
    
    # network styling
    fontsize = param.Selector(default = '10pt', objects = ['{}pt'.format(i) for i in range(31)])
    node_cmap = param.Selector(objects = hv.plotting.util.list_cmaps(), default='Spectral')
    node_clim = param.Tuple(default = (None, None))
    
    def __init__(self, 
                 nodes, 
                 edges,
                 index_col = 'GeneID',
                 source_col = 'GENE_ID_A',
                 target_col = 'GENE_ID_B',
                 label_col = 'geneSymbol',
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
        
        self.network_pane = pn.pane.HoloViews(sizing_mode = 'stretch_both', linked_axes=False, min_height=0, min_width=0)
        
        self.graph = DraggableGraph(
            index_col = self.index_col,
            source_col = self.source_col,
            target_col = self.target_col,
            label_col = self.label_col,
        )
        
        new_nodes, new_edges = self.update_nodes_edges(wait=True) # don't trigger param update yet
        
#         # set default style parameters if not passed in graph_opts
#         styles = list(zip(
#             ['Nodes', 'Graph', 'Nodes', 'Graph', 'Labels'], 
#             ['cmap', 'cmap', 'clim', 'clim', 'text_font_size'], 
#             ['node_cmap', 'node_cmap', 'node_clim', 'node_clim', 'fontize']
#         ))
        
#         for k1, k2, param_name in styles:
#             if not k1 in self.graph_opts:
#                 self.graph_opts[k1] = {}
#             if not k2 in self.graph_opts[k1]:
#                 self.graph_opts[k1].update({k2: getattr(self, param_name)})
                
        self.param.set_param(node_data=new_nodes, edge_data = new_edges, graph_opts = self.graph_opts) # triggers self.update_data
                 
    @param.depends('network_data.data', watch=True) # this method is the money maker!
    def view(self):
        network_graph = self.graph.view(self.network_data.data).opts(active_tools=['point_draw'])
        
        self.click_stream.source = self.graph.edge_graph
        self.network_pane.object = network_graph
        
        if self.selected_node == (None, None): # only on initialization
            self.selected_node = tuple(self.node_data.iloc[0,:][[self.index_col, self.label_col]])
        
        # for loading spinner control
        self.loading = False
        
    @param.depends('node_data', 'edge_data', 'graph_opts', 'layout', 'bundle_graph_edges', watch=True)
    def update_data(self):
        
        # for loading spinner control
        self.loading = True
        
        new_data = [
            self.node_data, # nodes
            self.edge_data, # edges
            [getattr(opts, k)(**self.graph_opts[k]) for k in self.graph_opts],
            self.layout,
            {'Yes': True, 'No': False}[self.bundle_graph_edges]            
        ]
        
        self.network_data.update(data=new_data) # triggers self.view
        
    @param.depends('click_stream.x', 'click_stream.y', watch=True)
    def get_clicked_node(self):
        x, y = (self.click_stream.x, self.click_stream.y)
        graph_data = pd.DataFrame(self.graph.stream.data)
        
        sel = graph_data.loc[((graph_data.x-x).abs()+(graph_data.y-y).abs()).idxmin(), :]
        idx = sel[self.index_col]
        label = sel[self.label_col]
        self.delta = (graph_data.x-x).abs()+(graph_data.y-y).abs().min()
        self.selected_node = (idx, label)
        
    @param.depends('max_nodes', watch=True)
    def update_nodes_edges(self, wait=False):
        
        # for loading spinner control
        self.loading = True
        
        new_nodes = self.nodes.iloc[:self.max_nodes,:].copy()
        new_edges = self.edges[self.edges[self.source_col].isin(new_nodes[self.index_col])&self.edges[self.target_col].isin(new_nodes[self.index_col])].copy()
        
        new_nodes['connectivity'] = pd.concat([new_edges.groupby(self.source_col).size(), new_edges.groupby(self.target_col).size()], axis=1).sum(axis=1).reindex(new_nodes[self.index_col]).fillna(0).values
        new_edges['edge_width'] = scale(new_edges['combined score'], 0.25, 5)
        
        self.node_clim = (new_nodes['connectivity'].min(), new_nodes['connectivity'].max()) # triggers self.update_node_clim
        self.sel_nodes = new_nodes
        
        if wait == True:
            return new_nodes, new_edges
        else:
            self.param.set_param(node_data = new_nodes, edge_data=new_edges, graph_opts = self.graph_opts) # triggers self.update_data
        
    @param.depends('fontsize', watch=True)
    def update_fontsize(self):
        if not 'Labels' in self.graph_opts:
            self.graph_opts['Labels'] = {}
        self.graph_opts['Labels'].update({'text_font_size': self.fontsize})
        self.param.set_param(graph_opts = self.graph_opts)
    
    @param.depends('node_cmap', watch=True)
    def update_node_cmap(self):
        if not 'Nodes' in self.graph_opts:
            self.graph_opts['Nodes'] = {}
        self.graph_opts['Nodes'].update({'cmap': self.node_cmap})
        if not 'Graph' in self.graph_opts:
            self.graph_opts['Graph'] = {}
        self.graph_opts['Graph'].update({'cmap': self.node_cmap})
        
        self.param.set_param(graph_opts = self.graph_opts)
        
    @param.depends('node_clim', watch=True)
    def update_node_clim(self): # Note: graph_opts must be set after call (e.g. after any change to self.node_clim)
        if not 'Nodes' in self.graph_opts:
            self.graph_opts['Nodes'] = {}
        self.graph_opts['Nodes'].update({'clim': self.node_clim})
        if not 'Graph' in self.graph_opts:
            self.graph_opts['Graph'] = {}
        self.graph_opts['Graph'].update({'clim': self.node_clim})