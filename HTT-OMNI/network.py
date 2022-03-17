import holoviews as hv
import param
import panel as pn

from draggable_graph import DraggableGraph

class Network(param.Parameterized):
    # will contain a list of [nodes, edges, graph_opts, layout_algorithm, bundle_graph_edges] to be displayed in the network
    network_data = hv.streams.Pipe(linked=False)
    
    # slider for # of nodes to show
    slider = param.Integer()
    
    def __init__(self, 
                 nodes, 
                 edges, 
                 graph_opts = [], 
                 layout_algorithm = 'circular', 
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
        self.network_data.update(data = [nodes, edges, graph_opts, layout_algorithm, bundle_graph_edges])
        
        self.slider = nodes.shape[0]
        self.param.slider.bounds = (1, nodes.shape[0])
        
    @param.depends('network_data.data', watch=True)
    def view(self):
        self.network_pane.object = self.graph.view(self.network_data.data)
        
    @param.depends('slider', watch=True)
    def update_nodes_edges(self):
        new_nodes = self.nodes.iloc[:self.slider,:]
        new_edges = self.edges[self.edges[self.source_col].isin(new_nodes[self.index_col])&self.edges[self.target_col].isin(new_nodes[self.index_col])].copy()
        
        new_data = [new_nodes, new_edges]+self.network_data.data[2:]
        self.network_data.update(data=new_data)