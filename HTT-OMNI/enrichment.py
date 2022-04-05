import param
import panel as pn
import holoviews as hv
import numpy as np
import pandas as pd
import requests
import seaborn as sns
from bokeh.models import HoverTool
from holoviews import opts, dim

class Enrichment(param.Parameterized):
    
    # enrichment query params
    results = param.DataFrame()
    selected_results = param.DataFrame()
    run_GO_analysis = param.Action(lambda x: x.param.trigger('run_GO_analysis'))
    GO_annot = param.Selector()
    
    # plot params
    GO_show = param.Integer(default = 10, bounds = (1, 50))
    GO_max_FDR = param.Number(0.05, bounds = (0, 1))
    GO_min_enrichment = param.Number(1, bounds = (1, 100))
    GO_plot_title = param.String()
    scatter_sizes = param.List([5, 10, 15, 20, 25, 30])
    scatter_bins = param.List([0, 10, 50, 100, 200, 500])
    
    def __init__(self, 
                 parent, 
                 annot_description_mapping, 
                 default_annot_query = None, 
                 index_col = 'geneID', 
                 ignore_terms = [], 
                 init_results = None,
                 **params
                ):
        super(Enrichment, self).__init__(**params)
        
        self.parent = parent
        self.index_col = index_col
        self.ignore_terms = ignore_terms
        self.mapping = annot_description_mapping
        
        # add annotation description mapping keys
        objects = np.sort(list(annot_description_mapping.keys()))
        self.param.GO_annot.objects = objects
        
        if default_annot_query is None:
            self.GO_annot = objects[0]
        else:
            self.GO_annot = default_annot_query
        
        self.plot_pane = pn.pane.HoloViews(sizing_mode='stretch_both', linked_axes=False, min_height=0, min_width=0)
        self.plot_pane.object = hv.DynamicMap(self.plot_GO_enrichment, streams = [self.param.selected_results, self.param.GO_min_enrichment, self.param.GO_max_FDR, self.param.GO_show])
        
        if init_results is None:
            self.param.trigger('run_GO_analysis')
        else:
            self.param.set_param(results=init_results)
        
        self.update_results()
            
    @param.depends('run_GO_analysis', watch=True)   
    def PantherGO_enrichment(self):
        
        # for loading spinner control
        self.loading = True
        
        url = 'http://pantherdb.org/services/oai/pantherdb/enrich/overrep'
        payload = {
            'geneInputList': ','.join(self.parent.sel_nodes[self.index_col].astype(str)), # comma-separated
            'organism': '9606',
            'annotDataSet': self.mapping[self.GO_annot]
        }

        r = requests.post(url, data=payload)
        results = r.json()['results']['result']
        results = pd.concat(list(map(pd.Series, results)), axis=1).T
        results = results[results['plus_minus']=='+']

        if len(results)>0:
            results = pd.concat([results, results['term'].apply(lambda x: pd.Series(x))], axis=1).drop('term', axis=1)
            results = results[~results['label'].isin(self.ignore_terms)]
        
        # for loading spinner control
        self.loading = False
        
        self.results = results.infer_objects().copy()
         
    def plot_GO_enrichment(self, selected_results, GO_min_enrichment, GO_max_FDR, GO_show):
        
        results = selected_results.copy()
        results['size'] = pd.cut(results['number_in_list'], bins=self.scatter_bins, labels = self.scatter_sizes[:-1]).astype(float).fillna(self.scatter_sizes[-1]).astype(int)
        results['num_out_of'] = results['number_in_list'].astype(str)+' (out of {}'.format(self.parent.sel_nodes[self.index_col].unique().shape[0])+')'


        cm = sns.blend_palette(['#ab4444', '#4489ab'], as_cmap=True)            
        clim = (results['fdr'].min(), GO_max_FDR)
        results = results.iloc[-GO_show:, :]

        hover = HoverTool(tooltips=[('Term', '@label'), 
                                ('Fold enrichment', '@fold_enrichment'), 
                                ('FDR', '@fdr'), 
                                ('# genes', '@num_out_of')]
                     )

        scatter = hv.Scatter(results, vdims=['fold_enrichment', 'size', 'fdr', 'num_out_of'], kdims = ['label'])
        spikes = hv.Spikes(results, kdims = ['label'], vdims = ['fold_enrichment', 'fdr'])

        scatter_opts = opts.Scatter(
            ylim = (None, results['fold_enrichment'].max()*1.1),
            size=dim('size'), 
            color='fdr', 
            cmap=cm, 
            responsive=True,
            ylabel = 'Fold enrichment',
            xlabel = '',
            tools = [hover],
            colorbar=True, 
            clim = clim, 
            colorbar_position = 'right', 
            clabel='False Discovery Rate',
            cnorm = 'log',
            invert_axes = True,
            shared_axes=False,
        )

        spike_opts = opts.Spikes(
            color='fdr', 
            line_width=3,
            cmap = cm, 
            cnorm = 'log', 
            clim=clim,
            invert_axes=True,
            toolbar='below',
            responsive=True,
            shared_axes=False,
        )
        
        return (spikes*scatter).opts(spike_opts, scatter_opts)

        
    @param.depends('results', 'GO_show', 'GO_min_enrichment', 'GO_max_FDR', watch=True)
    def update_results(self):
        n_tot = (self.results['fdr']<=self.GO_max_FDR).sum()
        selected_results = self.results[(self.results['fdr']<=self.GO_max_FDR)&(self.results['fold_enrichment']>=self.GO_min_enrichment)].sort_values(['fold_enrichment']).iloc[-self.GO_show:, :]

        if n_tot>0:
            self.GO_plot_title = 'Showing {} of {} significantly enriched {} terms'.format(selected_results.shape[0], n_tot, self.GO_annot,)
        else:
            self.GO_plot_title = 'No enriched terms'
        
        if selected_results.shape[0]>0:
            self.selected_results = selected_results
            self.plot_pane.visible = True
        else:
            self.plot_pane.visible = False
