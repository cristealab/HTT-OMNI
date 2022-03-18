import holoviews as hv
import param
import panel as pn

from draggable_graph import DraggableGraph

class Network(param.Parameterized):
    
    node_data = param.DataFrame()
    edge_data = param.DataFrame()
    graph_opts_data = param.Parameter()
    
    # data pipes to push data to DynamicMaps
    network_data = hv.streams.Pipe() #ultimately this will be a list of node_data, edge_data, graph_opts_data, layout, bundle_edge_graphs
    enrichment_data = hv.streams.Pipe()
    
    # slider for # of nodes to show
    max_nodes = param.Selector(objects = [10, 20, 50, 100, 200, 300, 400, 500], default=50)
    
    # graph layout algorithm
    layout = param.Selector(objects = ['kamada_kawai', 'circular', 'spring'], default='kamada_kawai')
    
    # edge bundling
    bundle_graph_edges = param.Selector(objects = ['Yes', 'No'], default='No')
    
    def __init__(self, 
                 nodes, 
                 edges, 
                 graph_opts = [], 
                 bundle_graph_edges = False, 
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
        self.network_data.update(data = [nodes, edges, graph_opts, self.layout, bundle_graph_edges])
        
        
    @param.depends('network_data.data', watch=True) # this method is the money maker!
    def view(self):
        self.network_pane.object = self.graph.view(self.network_data.data)
        
    @param.depends('max_nodes', watch=True)
    def update_nodes_edges(self):
        new_nodes = self.nodes.iloc[:self.max_nodes,:]
        new_edges = self.edges[self.edges[self.source_col].isin(new_nodes[self.index_col])&self.edges[self.target_col].isin(new_nodes[self.index_col])].copy()
        
        new_data = [new_nodes, new_edges]+self.network_data.data[2:]
        self.network_data.update(data=new_data)
        
    @param.depends('layout', watch=True)
    def update_layout(self):
        new_data = self.network_data.data[:3]+[self.layout, self.network_data.data[-1]]
        self.network_data.update(data=new_data)
        
    @param.depends('bundle_graph_edges', watch=True)
    def update_bundling(self):
        new_data = self.network_data.data[:4]+[{'Yes': True, 'No': False}[self.bundle_graph_edges]]
        self.network_data.update(data=new_data)
        