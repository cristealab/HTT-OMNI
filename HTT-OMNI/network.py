import holoviews as hv
import param
import panel as pn

from draggable_graph import DraggableGraph

class Network(param.Parameterized):
    
    node_data = param.DataFrame()
    edge_data = param.DataFrame() # list of [nodes, edges]
    graph_opts = param.Parameter({}) # dict of {k: v} where k is a valid opts.k method and v is a dict of options
    
    # data pipes to push data to DynamicMaps
    network_data = hv.streams.Pipe() #ultimately this will be a list of node_data, edge_data, graph_opts_data, layout, bundle_edge_graphs
    enrichment_data = hv.streams.Pipe()
    
    # slider for # of nodes to show
    max_nodes = param.Selector(objects = [10, 20, 50, 100, 200, 300, 400, 500], default=50)
    
    # graph layout algorithm
    layout = param.Selector(objects = ['kamada_kawai', 'circular', 'spring'], default='kamada_kawai')
    
    # edge bundling
    bundle_graph_edges = param.Selector(objects = ['Yes', 'No'], default='No')
    
    # network styling
    fontsize = param.Integer(default = 8, bounds = (1, 30))
    
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
        self.update_nodes_edges()
        
    @param.depends('network_data.data', watch=True) # this method is the money maker!
    def view(self):
        self.network_pane.object = self.graph.view(self.network_data.data).opts(active_tools=['point_draw'])
        
    @param.depends('node_data', 'edge_data', 'graph_opts', 'layout', 'bundle_graph_edges', watch=True)
    def update_data(self):
        new_data = [
            self.node_data, # nodes
            self.edge_data, # edges
            [getattr(opts, k)(**self.graph_opts[k]) for k in self.graph_opts],
            self.layout,
            {'Yes': True, 'No': False}[self.bundle_graph_edges]            
        ]
        
        self.network_data.update(data=new_data) # triggers self.view
        
    @param.depends('max_nodes', watch=True)
    def update_nodes_edges(self):
        new_nodes = self.nodes.iloc[:self.max_nodes,:]
        new_edges = self.edges[self.edges[self.source_col].isin(new_nodes[self.index_col])&self.edges[self.target_col].isin(new_nodes[self.index_col])].copy()
        
        self.param.set_param(node_data = new_nodes.copy(), edge_data=new_edges.copy())
        
    @param.depends('fontsize', watch=True)
    def update_fontsize(self):
        if not 'Labels' in self.graph_opts:
            self.graph_opts['Labels'] = {}
        self.graph_opts['Labels'].update({'text_font_size': '{}pt'.format(self.fontsize)})
        self.param.set_param(graph_opts = self.graph_opts)
        
        