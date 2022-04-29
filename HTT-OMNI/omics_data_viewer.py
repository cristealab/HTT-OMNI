import param
import pandas as pd
import holoviews as hv
from holoviews import opts
import numpy as np
import panel as pn

from network import Network

class OmicsDataViewer(param.Parameterized):
    data = param.DataFrame(precedence=-1) # subset of input data corresponding to selected ages/tissues
    selected_node_data = param.Series(precedence=-1) # subset of data corresponding to the selected node
    averaged_nodes_data = param.Series(precedence=-1) # subset of data corresponding to the network nodes (averaged)
    model_count = param.DataFrame(precedence=-1) # model count for all nodes in parent.nodes

    ages = param.ListSelector()
    tissues = param.ListSelector()
    
    parent = param.ClassSelector(Network, precedence=-1)
    net_node_ids = param.List()
    
    # plotting options
    plot_opts = param.Dict({'PROTEIN':[], 'RNA': [], 'SCRNA': [], 'SNRNA': []}, precedence=-1)
    type_protein = param.String('PROTEIN', constant=True, precedence=-1)
    type_RNA = param.String('RNA', constant=True, precedence=-1)
    type_snRNA = param.String('SNRNA', constant=True, precedence=-1)
    type_scRNA = param.String('SCRNA', constant=True, precedence=-1)
    type_model_obs = param.String('model_obs', constant=True, precedence=-1)
    title_network = param.String('All filtered nodes', constant=True, precedence=-1)
    
    def __init__(self, omics_data, dummy_leg, **params):
        super(OmicsDataViewer, self).__init__(**params)
        
        self.omics_data = omics_data
        self.dummy_leg = dummy_leg

        self.AS_types = ['PROTEIN', 'RNA']
        self.check_data()
        
        # widget mapping
        self.mapping = dict([
            ('tissues', {'type': pn.widgets.MultiChoice, 'name': 'Mouse AS tissues'}
            ), 
            ('ages', {'type': pn.widgets.MultiChoice, 'name': 'Mouse AS ages (months)'}
            ),
            ('age_label', {'type': pn.widgets.StaticText}
            ),
            ('tissue_label', {'type': pn.widgets.StaticText}
            )
        ])
        
        tissues_ = np.unique(self.omics_data[self.AS_types].columns.get_level_values('tissue')).tolist()
        ages_ = np.unique(self.omics_data[self.AS_types].columns.get_level_values('age')).tolist()
        
        self.param.tissues.objects = tissues_
        self.param.ages.objects = ages_
        
        self.param.set_param(tissues = tissues_, ages = ages_)
        self.update_data()
        
        pcts = np.array([0.05, 0.95])
        self.desc_lims = self.data.describe(percentiles=pcts).loc[['{:.0f}%'.format(pct) for pct in pcts*100], :]
        
        self.count_models()
        self.poly_overlays()
        
        # make plots
        pane_params = {'linked_axes': False, 'sizing_mode': 'stretch_both', 'min_height': 0, 'min_width': 0}
        
        self.selected_node_AS_protein_pane = pn.pane.HoloViews(**pane_params)
        stream1 = {'data': self.param.selected_node_data, 
                  'data_type': self.param.type_protein, 
                  'title': self.parent.param.selected_node}
        self.selected_node_AS_protein_pane.object = hv.DynamicMap(self.AS_plot, streams = stream1)
        
        self.selected_node_AS_RNA_pane = pn.pane.HoloViews(**pane_params)
        stream2 = {'data': self.param.selected_node_data, 
                  'data_type': self.param.type_RNA, 
                  'title': self.parent.param.selected_node}
        self.selected_node_AS_RNA_pane.object = hv.DynamicMap(self.AS_plot, streams = stream2)
        
        self.network_AS_protein_pane = pn.pane.HoloViews(**pane_params)
        stream3 = {'data': self.param.averaged_nodes_data, 
                  'data_type': self.param.type_protein, 
                  'title': self.param.title_network}
        self.network_AS_protein_pane.object = hv.DynamicMap(self.AS_plot, streams = stream3)
        
        self.network_AS_RNA_pane = pn.pane.HoloViews(**pane_params)
        stream4 = {'data': self.param.averaged_nodes_data, 
                  'data_type': self.param.type_RNA, 
                  'title': self.param.title_network}
        self.network_AS_RNA_pane.object = hv.DynamicMap(self.AS_plot, streams = stream4)
        
        self.selected_node_scRNA_pane = pn.pane.HoloViews(**pane_params)
        stream5 = {'data': self.param.selected_node_data, 
                  'data_type': self.param.type_scRNA, 
                  'title': self.parent.param.selected_node}
        self.selected_node_scRNA_pane.object = hv.DynamicMap(self.scRNA_plot, streams = stream5)
        
        self.network_scRNA_pane = pn.pane.HoloViews(**pane_params)
        stream6 = {'data': self.param.averaged_nodes_data, 
                  'data_type': self.param.type_scRNA, 
                  'title': self.param.title_network}
        self.network_scRNA_pane.object = hv.DynamicMap(self.scRNA_plot, streams = stream6)
        
        self.selected_node_snRNA_pane = pn.pane.HoloViews(**pane_params)
        stream7 = {'data': self.param.selected_node_data, 
                  'data_type': self.param.type_snRNA, 
                  'title': self.parent.param.selected_node}
        self.selected_node_snRNA_pane.object = hv.DynamicMap(self.snRNA_plot, streams = stream7)
        
        self.network_snRNA_pane = pn.pane.HoloViews(**pane_params)
        stream8 = {'data': self.param.averaged_nodes_data, 
                  'data_type': self.param.type_snRNA, 
                  'title': self.param.title_network}
        self.network_snRNA_pane.object = hv.DynamicMap(self.snRNA_plot, streams = stream8)
        
        self.selected_node_model_obs_pane = pn.pane.HoloViews(**pane_params)
        stream9 = {'ids': self.parent.param.selected_node, 
                  'data_': self.param.model_count,
                  'data_type': self.param.type_model_obs, 
                  'title': self.param.title_network} # cannot be self.parent.param.selected_node
        self.selected_node_model_obs_pane.object = hv.DynamicMap(self.models_plot, streams = stream9)
        
        self.network_model_obs_pane = pn.pane.HoloViews(**pane_params)
        stream10 = {'ids': self.param.net_node_ids, 
                  'data_': self.param.model_count,
                  'data_type': self.param.type_model_obs, 
                  'title': self.param.title_network}
        self.network_model_obs_pane.object = hv.DynamicMap(self.models_plot, streams = stream10)
        
    def check_data(self):
        if not self.omics_data.columns.names == ['type', 'tissue', 'Q-length', 'age', 'Tissue/Age']:
            raise ValueError('Column index names do not match "[type, tissue, Q-length, age, Tissue/Age]"')
            
        if not self.omics_data.index.names == ['geneID', 'geneSymbol']:
            raise ValueError('Index names must match [geneID, geneSymbol]')
        
        if self.omics_data.columns.get_level_values('Q-length').dtype!=np.int64:
            raise ValueError('Q-length column values must be int64 type, got type {}'.format(self.omics_data.columns.get_level_values('Q-length').dtype))
        
        if self.omics_data.columns.get_level_values('age').dtype!=np.int64:
            raise ValueError('Age column values must be int64 type, got type {}'.format(self.omics_data.columns.get_level_values('age').dtype))
        
        if self.omics_data.index.duplicated().any():
            raise ValueError('Index contains duplicate entries')
        
        if any([self.omics_data.columns.get_level_values(l).isnull().any() for l in self.omics_data.columns.names]):
            raise ValueError('Column index contains null values')
            
        if any([self.omics_data.index.get_level_values(l).isnull().any() for l in self.omics_data.index.names]):
            raise ValueError('Index contains null values')
    
    @param.depends('parent.parent.nodes', watch = True)
    def count_models(self):
        self.model_count = self.parent.parent.nodes.groupby(self.parent.parent.groupby_PPI_cols+['model']).size().groupby([self.parent.parent.index_col,'model']).size().unstack('model').fillna(0)

    @param.depends('ages', 'tissues', watch=True)
    def update_data(self):
        
        cols = self.omics_data.columns
        data_ = self.omics_data[cols[(cols.get_level_values('tissue').isin(self.tissues)&cols.get_level_values('age').isin(self.ages))|(~cols.get_level_values('type').isin(self.AS_types))]]
        
        self.data = data_
    
    @param.depends('parent.selected_node', 'data', watch=True)
    def update_selected_node_data(self):
        
        data_ = self.data.reindex([self.parent.selected_node]).mean()
        data_.name = 'value'
        
        self.selected_node_data = data_
    
    @param.depends('parent.sel_nodes', 'data', watch=True)
    def update_averaged_nodes_data(self):
        idx = list(map(tuple, self.parent.sel_nodes[[self.parent.index_col, self.parent.label_col]].values))
        data_ = self.data.reindex(idx).mean()
        data_.name = 'value'
        self.param.set_param(averaged_nodes_data = data_, net_node_ids = self.parent.sel_nodes[self.parent.index_col].values.tolist())
    
    def AS_plot(self, data, data_type, title):
        data_ = data[data_type].reset_index()
        ds = hv.Dataset(data_)
        
        if type(title) == tuple:
            title_ = title[1]
            ylabel = 'log2FC abundance'
        else:
            title_ = title
            ylabel = 'log2FC abundance (mean of all nodes)'

        lines = ds.to(hv.Curve, 'Q-length', ['value', 'tissue','age'], 'Tissue/Age').overlay('Tissue/Age')
        scatter = ds.to(hv.Scatter, 'Q-length', ['value', 'tissue','age'], 'Tissue/Age').overlay('Tissue/Age')

        return (self.polys[data_type]*self.dummy_leg*((lines*scatter).opts(self.plot_opts[data_type]))).opts(title=title_, ylabel=ylabel, framewise=True)
    
    def scRNA_plot(self, data, data_type, title):
        data_ = data[data_type].reset_index()[['tissue', 'value']]
        data_.columns = ['cell type', 'fractional expression']
        data_['alpha'] = np.where(data_['fractional expression']>=0.95, 1, 0.15)
        ds = hv.Dataset(data_)

        if type(title) == tuple:
            title_ = title[1]
            ylabel = 'fractional expression'
        else:
            title_ = title
            ylabel = 'fractional expression (mean of all nodes)'

        return hv.Bars(ds, 'cell type', ['fractional expression', 'alpha']).opts(self.plot_opts[data_type]).opts(ylabel=ylabel, title=title_, framewise=True)
    
    def snRNA_plot(self, data, data_type, title):
        data_ = data[data_type].reset_index()[['tissue', 'value']]
        data_.columns = ['cell type', 'log2FC abundance (Q175/WT)']
        sn_lims = self.desc_lims[data_type]
        sn_lims = sn_lims.T.reset_index([i for i in sn_lims.columns.names if not i=='tissue'], drop=True).T

        sig = (data_['log2FC abundance (Q175/WT)']<=sn_lims.iloc[0,:][data_['cell type']].values)|(data_['log2FC abundance (Q175/WT)']>=sn_lims.iloc[1,:][data_['cell type']].values)
        data_['alpha'] = np.where(sig, 1, 0.15)
        ds = hv.Dataset(data_)

        if type(title) == tuple:
            title_ = title[1]
            ylabel = 'log2FC abundance (Q175/WT)'
        else:
            title_ = title
            ylabel = 'log2FC abundance\n(Q175/WT; mean of all nodes)'

        return hv.Bars(ds, 'cell type', ['log2FC abundance (Q175/WT)', 'alpha']).opts(self.plot_opts[data_type]).opts(ylabel=ylabel, title=title_)
    
    def models_plot(self, ids, data_, data_type, title):
        if isinstance(ids, tuple):
            data = data_.loc[ids[0], :].reset_index()
            title_ = ids[1]
            ylabel = '# PPI observations'
        elif isinstance(ids, list):
            data = data_.reindex(ids).sum().reset_index()
            title_ = title
            ylabel = '# PPI observations (sum of all nodes)'

        data.columns = ['model', '# PPI observations']

        return hv.Bars(data).opts(self.plot_opts[data_type]).opts(title=title_, ylabel=ylabel)
    
    def poly_overlays(self):
        lims = self.desc_lims[self.AS_types].stack(['type', 'Q-length']).mean(axis=1).unstack(['type', 'Q-length'])
        poly_opts = opts.Polygons(
            color='black', 
            line_color=None, 
            fill_alpha=0.05,
            shared_axes=False,
            show_legend=True,
            xlabel = 'Q-length'
        )

        self.polys = {}
        for t in self.AS_types:
            xs = lims[t].iloc[0].index.values.tolist()+lims[t].iloc[1].index.values.tolist()[::-1]
            ys = lims[t].iloc[0].values.tolist()+lims[t].iloc[1].values.tolist()[::-1]
            poly = hv.Polygons([{'x': xs, 'y': ys}])
            self.polys[t] = poly.opts(poly_opts)