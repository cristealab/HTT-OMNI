import pandas as pd
import numpy as np
import holoviews as hv
import networkx as nx
import param
from holoviews import opts, dim
from holoviews.operation.datashader import bundle_graph
import panel as pn

hv.extension('bokeh')
pn.extension(sizing_mode='stretch_width', min_height=0, min_width=0)

class DraggableGraph(param.Parameterized):
    
    nodes = param.DataFrame(precedence=-1)
    edges = param.DataFrame(precedence=-1)
    
    node_data = param.DataFrame()
    
    stream = hv.streams.PointDraw(add=False)
    
    def __init__(self, 
                 nodes, 
                 edges, 
                 index_col = 'geneID', 
                 source_col = 'GENE_ID_A', 
                 target_col = 'GENE_ID_B', 
                 label_col = 'geneSymbol',
                 graph_opts = [],
                 bundle_graph_edges = False,
                 layout_algorithm = 'circular',
                 **params
                ):
        
        super(DraggableGraph, self).__init__(**params)
        
        self.index_col = index_col
        self.source_col = source_col
        self.target_col = target_col
        self.label_col = label_col
        self.bundle_graph_edges = bundle_graph_edges
        
        # make sure that index, source, and target are the same dtype
        if np.unique([nodes.dtypes[index_col], edges.dtypes[source_col], edges.dtypes[target_col]]).shape[0]>1:
            raise TypeError('Index, source, and target columns must have the same dtype')
            
        # make sure that all source and target values are present in the index column
        # (i.e., that all edges have nodes)
        if not (edges[self.source_col].isin(nodes[self.index_col]).all()&edges[self.target_col].isin(nodes[self.index_col]).all()):
            raise TypeError('Values in source and/or target columns not present in node index column')
        
        self.G = self.make_graph(nodes, edges)
        
        self.nodes = nodes
        self.edges = edges
        
        init_layout = pd.DataFrame(getattr(nx, '{}_layout'.format(layout_algorithm))(self.G), index=['x', 'y']).T
        init_layout.index.name = index_col
        
        self.node_data = pd.concat([self.nodes.set_index(index_col), init_layout], axis=1).reset_index()
        
        self.node_graph = hv.DynamicMap(self.view_nodes)
        self.stream.source = self.node_graph
        
        self.edge_graph = hv.DynamicMap(self.view_edges, streams = [self.stream])
        self.labels = hv.DynamicMap(self.view_labels, streams = [self.stream])
        
        overlay = (self.edge_graph*self.node_graph*self.labels)        
        self.overlay = pn.pane.HoloViews(overlay.opts(*graph_opts), 
                                         sizing_mode='stretch_both', 
                                         height=400, 
                                         width=400
                                         )
    
    def make_graph(self, nodes, edges):
        
        G = nx.Graph()
        for idx, data in nodes.sort_index().iterrows():
            G.add_node(data[self.index_col], **data.to_dict())
        
        for idx, data in edges.sort_index().iterrows():
            G.add_edge(data[self.source_col], data[self.target_col], **data.to_dict())
         
        return G
        
    def view_nodes(self):
        g = hv.Graph.from_networkx(self.G, dict(zip(self.node_data[self.index_col], self.node_data[['x', 'y']].values)))
                
        return g.nodes
        
    def view_edges(self, data):
        if data is None:
            data = self.node_data
        else:
            data = pd.DataFrame(data)
        
        g = hv.Graph.from_networkx(self.G, dict(zip(data[self.index_col], data[['x', 'y']].values)))
        
        if self.bundle_graph_edges == True:
            g = bundle_graph(g)
        
        return g
    
    def view_labels(self, data):
        if data is None:
            data = self.node_data
        else:
            data = pd.DataFrame(data)
        
        return hv.Labels(data, ['x', 'y'], self.label_col)
        
    @param.depends('stream.data', watch=True)
    def _update_dataframe(self):
        self.node_data = pd.DataFrame(self.stream.data)