from turtle import back
import holoviews as hv
import panel as pn
import gc

from data_filter import DataFilter
from network import Network
from enrichment import Enrichment
from omics_data_viewer import OmicsDataViewer
from app import App

# from memory_profiler import profile

hv.extension('bokeh')

# @profile
def user_instance():
    
    gc.collect()

    data_filter = DataFilter(**{k:pn.state.cache[k] for k in ['nodes', 'edges', 'filters', 'index_col', 'gene_symbol_col', 'filter_aliases', 'groupby_PPI_cols']})

    network = Network(parent = data_filter, 
                      graph_opts = pn.state.cache['graph_opts'].copy(), 
                      **{k:pn.state.cache[k] for k in ['nodes', 'edges', 'index_col', 'source_col', 'target_col', 'label_col', 'fontsize', 'node_cmap', 'user_tooltips']})

    enrichment = Enrichment(parent = network, **{k:pn.state.cache[k] for k in ['annot_description_mapping', 'index_col', 'background_geneIDs']})  

    omics_viewer = OmicsDataViewer(parent = network, **{k:pn.state.cache[k] for k in ['omics_data', 'dummy_leg', 'plot_opts']})

    app = App(
        enrichment = enrichment, 
        data_filter=data_filter, 
        network=network, 
        omics_viewer=omics_viewer,
    )

    return app.view()

if __name__ == '__main__':
    from config_setup import setup
    setup() # run setup file to read variables into pn.state.cache
    pn.serve(user_instance, show = True)
else:
    user_instance()