import holoviews as hv
import pandas as pd
import numpy as np

# legends needed by the app:
### colorbar for the network
### node sizes for the network
### scatter sizes for the GO plot
### Tissue/age for the omics plots


def make_colorbar(cmap, clim, orientation = 'horizontal', position = 'top', **kwargs):
    hm = hv.HeatMap([(0, 0, clim[0]), (0, 1, clim[1])])

    hm_opts = {
        'colorbar': True, 
        'colorbar_opts': kwargs,
        'alpha': 0,
        'responsive': True, 
        'colorbar_position': position,
        'cmap': cmap,
        'xaxis': None,
        'yaxis': None,
        'show_frame': True,
        'bgcolor': (0, 0, 0, 0),
        'frame_height': 1,
    }

    return hm.opts(**hm_opts)

def nodes_colorbar(cmap, clim):
    return make_colorbar(cmap, 
                         clim, 
#                           height = 18, 
                          #width = 'auto', 
                          title = 'node color = connectivity', 
                          background_fill_alpha=0, 
                          padding=0, 
                          )

def tissue_age_leg(all_tissues, all_ages, x, graph_opts = []):
    '''
    dummy_leg_opts = [
        opts.Curve(
            line_color = dim('index').categorize(tissue_colors),
            line_dash=dim('index').categorize(age_dashes),
            framewise=True,    
        ),
        opts.NdOverlay(
            legend_position='right', 
            show_legend=True,
            legend_opts = {'title':'Tissue/Age (mo)'},
            framewise=True,
            show_frame = False,
            hooks = [save_hook]
        )
    ]

    dummy_leg = tissue_age_leg(all_tissues, all_ages, 150, graph_opts = dummy_leg_opts)

    '''

    tissue_age = all_tissues.tolist()+all_ages.tolist()
    
    dummy = pd.DataFrame(np.zeros(shape=(len(tissue_age), 2)), index=tissue_age, columns=[x, x])
    dummy_ds = hv.Dataset(dummy.reset_index().melt('index'))
    dummy_ds.columns = ['tissue_age', 'variable', 'value']

    return dummy_ds.to(hv.Curve, 'variable', ['value'], 'index').overlay('index').opts(*graph_opts)
