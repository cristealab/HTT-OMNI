import holoviews as hv

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
