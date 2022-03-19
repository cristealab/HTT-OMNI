import pandas as pd
import numpy as np
import holoviews as hv
import networkx as nx
import param
from holoviews.operation.datashader import bundle_graph

class DraggableGraph(param.Parameterized):
    
    # keeps track of previous nodes, edges, and layout to maintain node positions when changing aesthetic properties
    current_nodes = param.DataFrame(precedence=-1)
    current_edges = param.DataFrame(precedence=-1)
    current_layout = param.String()
    current_stream_data = param.Parameter()
    
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
        
    def view_nodes(self, positions):
        g = hv.Graph.from_networkx(self.G, dict(zip(positions[self.index_col], positions[['x', 'y']].values))).nodes
        self.stream.update(data=g.columns())

        return g
        
    def view_edges(self, data):

        data_ = pd.DataFrame(data)
        
        g = hv.Graph.from_networkx(self.G, dict(zip(data_[self.index_col], data_[['x', 'y']].values)))
        
        if self.bundle_graph_edges == True:
            g = bundle_graph(g)
        
        return g
    
    def view_labels(self, data):
        
        data_ = pd.DataFrame(data)
        
        return hv.Labels(data_, ['x', 'y'], self.label_col)
    
    def view(self, data):
        if len(data)!=5:
            raise ValueError('Data does not have the right number of items (nodes, edges, graph_opts, layout_algorithm, bundle_graph_edges)')
        
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
        
        if (self.current_nodes is None) and (self.current_edges is None):
            new_layout = True
        elif self.current_layout!=layout_algorithm:
            new_layout = True
        elif self.current_nodes.index.isin(nodes.index).all() and nodes.index.isin(self.current_nodes.index).all():
            new_layout = False
        else:
            new_layout = True
            
        if new_layout == True:
            init_layout = pd.DataFrame(getattr(nx, '{}_layout'.format(layout_algorithm))(self.G), index=['x', 'y']).T
            init_layout.index.name = self.index_col
            positions = pd.concat([nodes.set_index(self.index_col), init_layout], axis=1).reset_index()
        else:
            positions = pd.DataFrame(self.current_stream_data)
            
        self.node_graph = self.view_nodes(positions)
        self.stream.source = self.node_graph
        
        self.edge_graph = hv.DynamicMap(self.view_edges, streams = [self.stream])
        self.labels = hv.DynamicMap(self.view_labels, streams = [self.stream])
        
        self.current_nodes = nodes
        self.current_edges = edges
        self.current_layout = layout_algorithm
        
        return (self.edge_graph*self.node_graph*self.labels).opts(*graph_opts)        
    
    @param.depends('stream.data', watch=True)
    def _update(self):
        self.node_data = pd.DataFrame(self.stream.data)
        self.current_stream_data = self.stream.data
        