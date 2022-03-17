import pandas as pd
import numpy as np
import holoviews as hv
import networkx as nx
import param
from holoviews.operation.datashader import bundle_graph

class DraggableGraph(param.Parameterized):
    
    nodes = param.DataFrame(precedence=-1)
    edges = param.DataFrame(precedence=-1)
    
    node_data = param.DataFrame()
    stream = hv.streams.PointDraw(add=False)
    
    def __init__(self, 
                 index_col = 'geneID', 
                 source_col = 'GENE_ID_A', 
                 target_col = 'GENE_ID_B', 
                 label_col = 'geneSymbol',
                 **params
                ):
        
        super(DraggableGraph, self).__init__(**params)
        
        # holoviews throws an error when bundling edges if source and target == 'source' or 'target'
        # this is likely a bug (perhaps report later?)
        if source_col == 'source':
            raise ValueError('source_col cannot be "source" (passed "{}")'.format(source_col))
        elif target_col=='target':
            raise ValueError('target_col cannot be "target" (passed "{}")'.format(target_col))
        
        self.index_col = index_col
        self.source_col = source_col
        self.target_col = target_col
        self.label_col = label_col
            
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
            data_ = self.node_data
        else:
            data_ = pd.DataFrame(data)
        
        g = hv.Graph.from_networkx(self.G, dict(zip(data_[self.index_col], data_[['x', 'y']].values)))
        
        if self.bundle_graph_edges == True:
            g = bundle_graph(g)
        
        return g
    
    def view_labels(self, data):
        if data is None:
            data_ = self.node_data
        else:
            data_ = pd.DataFrame(data)
        
        return hv.Labels(data_, ['x', 'y'], self.label_col)
    
    def view(self, data):
        if len(data)!=5:
            raise ValueError('Data does not have the right number of items (nodes, edges, graph_opts, layout_algorithm, bundle_graph_edges)')
            
        # make a new stream (for some reason we need a new stream for new data...)
        self.stream = hv.streams.PointDraw(add=False)
        
        # unpack data = [nodes, edges, graph_opts, layout_algorithm, bundle_graph]
        nodes, edges, graph_opts, layout_algorithm, bundle_graph_edges = data
        
        # make sure that index, source, and target are the same dtype
        if np.unique([nodes.dtypes[self.index_col], edges.dtypes[self.source_col], edges.dtypes[self.target_col]]).shape[0]>1:
            raise TypeError('Index, source, and target columns must have the same dtype')
            
        # make sure that all source and target values are present in the index column
        # (i.e., that all edges have nodes)
        if not (edges[self.source_col].isin(nodes[self.index_col]).all()&edges[self.target_col].isin(nodes[self.index_col]).all()):
            raise TypeError('Values in source and/or target columns not present in node index column')
        
        self.G = self.make_graph(nodes, edges)
        self.bundle_graph_edges = bundle_graph_edges
        
        init_layout = pd.DataFrame(getattr(nx, '{}_layout'.format(layout_algorithm))(self.G), index=['x', 'y']).T
        init_layout.index.name = self.index_col
        
        self.node_data = pd.concat([nodes.set_index(self.index_col), init_layout], axis=1).reset_index()
        
        node_graph = self.view_nodes()
        self.stream.source = node_graph
        
        edge_graph = hv.DynamicMap(self.view_edges, streams = [self.stream])
        labels = hv.DynamicMap(self.view_labels, streams = [self.stream])
        
        return (edge_graph*node_graph*labels).opts(*graph_opts)        
    
    @param.depends('stream.data', watch=True)
    def _update_dataframe(self):
        self.node_data = pd.DataFrame(self.stream.data)