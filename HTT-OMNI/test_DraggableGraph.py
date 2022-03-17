import pandas as pd
import numpy as np
import itertools
import holoviews as hv
from holoviews import opts, dim
import panel as pn

from draggable_graph import DraggableGraph

hv.extension('bokeh')
pn.extension(sizing_mode='stretch_width', min_height=0, min_width=0)
                
def make_nodes_edges(n, links=None):
    if links is None:
        links = n
        
    nodes = pd.DataFrame(list(zip(np.random.uniform(low=20, high=90, size=n), [str(i) for i in range(1, n+1)])), columns = ['node_score', 'label'])
    e = np.array(list(itertools.combinations(nodes['label'][:links], 2)))
    edges = pd.DataFrame(np.vstack([e.T, np.random.uniform(0, 256, len(e))]).T, columns = ['source_node', 'target_node', 'edge_score'])
    edges['edge_score'] = edges['edge_score'].astype(float)
    
    return nodes, edges

test_nodes, test_edges = make_nodes_edges(10)

node_color = 'node_score'
node_size = 'node_score'
node_cmap = 'Spectral'

graph_opts = [
    opts.Nodes(
        active_tools=['point_draw',], 
        color=node_color,
        cmap = node_cmap,
        size = node_size,
    ),
    opts.Graph(
        edge_color = dim('edge_score'),
        edge_cmap = 'Greys',
        node_color = node_color,
        cmap = node_cmap,
        node_size = node_size,
        tools = []
    ),
    opts.Overlay(
        xlim=(-1.3, 1.3),
        ylim=(-1.3, 1.3),
        responsive=True,
        xaxis=None, 
        yaxis=None
    )
]

graph = DraggableGraph(
    index_col = 'label',
    source_col = 'source_node',
    target_col = 'target_node',
    label_col = 'label',
)

pane = pn.pane.HoloViews(
    graph.view([test_nodes.copy(), test_edges.copy(), graph_opts, 'circular', False]), 
    sizing_mode='stretch_both', 
    linked_axes=False, 
    min_height=0, 
    min_width=0
)

react = pn.template.ReactTemplate()
react.main[:4,:6] = pn.Column(pane)
react.show()