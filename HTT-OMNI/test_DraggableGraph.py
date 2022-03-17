import pandas as pd
import numpy as np
import itertools
import holoviews as hv
from holoviews import opts, dim
import panel as pn

from draggable_graph import DraggableGraph

hv.extension('bokeh')
pn.extension(sizing_mode='stretch_width', min_height=0, min_width=0)
                
test_nodes = pd.DataFrame([[20., 'a'],
                           [50., 'b'],
                           [30., 'c'],
                           [50., 'd'],
                           [75., 'e'],
                           [88., 'f']], columns = ['node_score', 'label'])

e = np.array(list(itertools.combinations(test_nodes['label'], 2)))
test_edges = pd.DataFrame(np.vstack([e.T, np.random.uniform(0, 256, len(e))]).T, columns = ['source', 'target', 'edge_score'])
test_edges['edge_score'] = test_edges['edge_score'].astype(float)

graph_opts = [
    opts.Nodes(
        active_tools=['point_draw',], 
        color='white',
        line_color='black',
        size = 'node_score',
    ),
    opts.Graph(
        edge_color = dim('edge_score'),
        node_size = 0, # hide nodes from this representation
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
    nodes=test_nodes, 
    edges=test_edges,
    index_col = 'label',
    source_col = 'source',
    target_col = 'target',
    label_col = 'label',
    graph_opts = graph_opts,
)
react = pn.template.ReactTemplate()
react.main[:4,:4] = pn.Column(graph.overlay)
react.show()