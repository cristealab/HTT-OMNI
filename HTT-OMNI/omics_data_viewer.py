import param
import pandas as pd
import holoviews as hv
import numpy as np
import panel as pn

from network import Network

class OmicsDataViewer(param.Parameterized):
    data = param.DataFrame() # subset of input data corresponding to selected ages/tissues
    selected_node_data = param.Series() # subset of data corresponding to the selected node
    averaged_nodes_data = param.Series() # subset of data corresponding to the network nodes (averaged)
    
    ages = param.ListSelector()
    tissues = param.ListSelector()
    
    parent = param.ClassSelector(Network)
    
    # plotting options
    plot_opts = param.Dict({'Protein/RNA':[], 'scRNA': [], 'snRNA': []})
    type_protein = param.String('PROTEIN', constant=True)
    type_RNA = param.String('RNA', constant=True)
    title_network = param.String('All filtered nodes', constant=True)
    
    def __init__(self, omics_data, **params):
        super(OmicsDataViewer, self).__init__(**params)
        
        self.omics_data = omics_data
        self.AS_types = ['PROTEIN', 'RNA']
        self.check_data()
        
        tissues_ = np.unique(self.omics_data[self.AS_types].columns.get_level_values('tissue')).tolist()
        ages_ = np.unique(self.omics_data[self.AS_types].columns.get_level_values('age')).tolist()
        
        self.param.tissues.objects = tissues_
        self.param.ages.objects = ages_
        
        self.param.set_param(tissues = tissues_, ages = ages_)
        self.update_data()
        
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
        self.averaged_nodes_data = data_
    
    def AS_plot(self, data, data_type, title):
        data_ = data[data_type].reset_index()
        ds = hv.Dataset(data_)
        
        if type(title) == tuple:
            title_ = title[1]
        else:
            title_ = title

        lines = ds.to(hv.Curve, 'Q-length', ['value', 'tissue','age'], 'Tissue/Age').overlay('Tissue/Age')
        scatter = ds.to(hv.Scatter, 'Q-length', ['value', 'tissue','age'], 'Tissue/Age').overlay('Tissue/Age')

        return (lines*scatter).opts(self.plot_opts[data_type]).opts(title=title_, framewise=True)