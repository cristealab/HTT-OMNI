import param
import panel as pn
import pandas as pd
import numpy as np

class DataFilter(param.Parameterized):
    filters = param.List(precedence=-1)
    
    filtered_nodes = param.DataFrame(precedence=-1)
    queried_nodes = param.DataFrame(precedence=-1)
    sel_nodes = param.DataFrame(precedence=-1)
    sel_edges = param.DataFrame(precedence=-1)
    show_nodes = param.DataFrame(precedence=-1)
    show_edges = param.DataFrame(precedence=-1)
    
    # node params
    node_query = param.String(default='')
    max_nodes = param.Selector(objects = [10, 20, 50, 100, 200, 300, 400, 500, 1000, 4000], default=50)
    node_display_priority = param.Selector(objects = ["# PPI observations (all)", "# PPI observations (filtered)"], default = '# PPI observations (all)')
    vis_unconnected = param.Selector(objects = ['Hide', 'Show'], default='Show')

    # edge params
    STRINGdb_score = param.Number(0.4, bounds=(0, 1))
    
    def __init__(self, 
                 nodes, 
                 edges, 
                 index_col = 'geneID',
                 gene_symbol_col = 'geneSymbol',
                 source_col = 'GENE_ID_A',
                 target_col = 'GENE_ID_B',
                 edge_score_col = 'combined_score',
                 filter_aliases = None,
                 groupby_PPI_cols = ['geneID', 'studyID'],
                 **params):
        
        super(DataFilter, self).__init__(**params)
        
        self.nodes = nodes
        self.edges = edges
        self.index_col = index_col
        self.gene_symbol_col = gene_symbol_col
        self.groupby_PPI_cols = groupby_PPI_cols
        self.source_col = source_col
        self.target_col = target_col
        self.edge_score_col = edge_score_col
        
        if filter_aliases is None:
            filter_aliases = {k: k for k in self.filters}
        
        self.check_data()
        
        self.PPI_sum = self.compute_PPI_sum(nodes)
        self.sym_to_index = self.nodes[[self.index_col, self.gene_symbol_col]].drop_duplicates()
        self.one_hot = self.encode_one_hot(nodes)
        self.annotations = self.one_hot.groupby(level=0, axis=1).apply(self.one_hot_to_str).reset_index(self.gene_symbol_col)
        self.annotations['PPI_SUM_TOTAL'] = self.PPI_sum.reindex(self.annotations.index)
        
        options_ = self.nodes[self.filters].apply(lambda x: np.unique(x.dropna()).tolist()).to_dict()
        self.options_ = options_
        self.options_map = pd.concat({f: self.nodes.groupby([self.index_col, f]).size().unstack().notnull().apply(lambda x: '{} ({})'.format(x.name, x.sum())) for f in self.filters})
        self.options_map_r = pd.Series({(a, c): b for (a, b), c in self.options_map.iteritems()})
        
        for opt in self.options_:
            self.param._add_parameter(opt, param.ListSelector(default = [], objects = self.options_map[opt].values.tolist()))
            self.param._add_parameter(opt+'_AND_OR_NOT', param.Selector(default = 'OR', objects = ['AND', 'OR', 'NOT']))
            
        self.param.watch(self.filter_nodes, self.filters+[opt+'_AND_OR_NOT' for opt in self.options_])
        
        # widget mapping
        default = [(k, {'type': pn.widgets.MultiChoice, 'solid': False, 'placeholder': 'SHOW ALL', 'name': filter_aliases[k]}) if len(options_[k])<1000 else (k, {'type': pn.widgets.MultiSelect, 'size':10}) for k in options_]
        default_AND_OR_NOT = [(k+'_AND_OR_NOT', {'type': pn.widgets.RadioButtonGroup}) for k in options_]
        other = [
            ('node_query', {'type': pn.widgets.TextAreaInput, 
                                 'placeholder': 'Type gene IDs or gene symbols of interest (one per line) e.g.,\n123\n456\n789', 
                                 'max_length': 10**10, 
                                 'min_height': 100, 
                                 'sizing_mode': 'stretch_both'}
            ),
            ('max_nodes', {'type': pn.widgets.DiscreteSlider, 
                                'throttled': True}
            ),
            ('vis_unconnected', {'type': pn.widgets.RadioButtonGroup}
            ),
            ('STRINGdb_score', {'throttled': True}
            ),
        ]
        
        self.mapping = dict(other+default+default_AND_OR_NOT)

        self.filter_nodes(1)# triggers self.apply_query, self.update_sel_nodes

    def check_data(self):
        if not np.isin(self.filters, self.nodes.columns).all():
            missing = np.array(self.filters)[~np.isin(self.filters, self.nodes.columns)]
            raise KeyError('Not all filters were found as columns in nodes (missing: {})'.format(missing))
            
        if (self.nodes.groupby([self.index_col, self.gene_symbol_col]).size().groupby(self.index_col).size()>1).any():
            temp = self.nodes.groupby([self.index_col, self.gene_symbol_col]).size().groupby(self.index_col).size()>1
            offending_genes = temp.index[temp].values.tolist()            
            raise ValueError('Some gene IDs map to more than one gene symbol (offending gene IDs = {})'.format(offending_genes))
        
        if not self.index_col in self.groupby_PPI_cols:
            raise KeyError('"index_col" ({}) must be present in "groupby_PPI_cols" ({})'.format(self.index_col, self.groupby_PPI_cols))
            
    def encode_one_hot(self, df):
        one_hot = pd.concat({f: df.groupby([self.index_col, self.gene_symbol_col, f]).size().unstack() for f in self.filters}, axis=1).fillna(0).where(lambda x: x==0, 1).astype(bool)
        one_hot = one_hot[one_hot.columns[one_hot.columns.get_level_values(1)!='Not reported']].copy()
        
        return one_hot
    
    def one_hot_to_str(self, df):
        arr_str = np.where(df, df.columns.get_level_values(1), 'EMPTY')+', '

        return pd.Series(arr_str.sum(axis=1), index = df.index).str.replace('EMPTY, ', '').str.strip(', ')
    
    def compute_PPI_sum(self, df):
        return df.groupby(self.groupby_PPI_cols).size().groupby(self.index_col).size()
    
    def filter_nodes(self, event):
        activated_filters = pd.Series([f if getattr(self, f)!=[] else np.nan for f in self.filters], index = self.filters)

        if not activated_filters.isnull().all():
            current_filters = pd.Series([self.options_map_r[f][getattr(self, f)].values.tolist() if getattr(self, f)!=[] else np.nan for f in self.filters], index = self.filters)

            # filter nodes based on activated filter selection and AND/OR toggles
            temp = self.nodes.groupby([self.index_col]+activated_filters.dropna().values.tolist()).apply(lambda x: ' '.join(x.index.astype(str)))

            for f in activated_filters.dropna():
                bool_ = getattr(self, f+'_AND_OR_NOT')

                if bool_ in ['AND', 'OR']:
                    us = temp.unstack(f).T.reindex(current_filters[f]).T

                    if bool_ == 'AND':
                        temp = us[us.notnull().all(axis=1)].stack(f)
                    elif bool_ == 'OR':
                        temp = us[us.notnull().any(axis=1)].stack(f)

                elif bool_ == 'NOT':
                    idx = self.one_hot.index.get_level_values(self.index_col)[~(self.one_hot[f][current_filters[f]]==1).any(axis=1)]
                    temp = temp[temp.index.get_level_values(self.index_col).isin(idx)]

            if temp.shape[0]>0:
                filtered_nodes = self.nodes.loc[temp.str.split(' ', expand=True).stack().astype(int).values, :]
            else:
                filtered_nodes = self.nodes.iloc[[], :]
        
        else:
            filtered_nodes = self.nodes
        
        self.filtered_nodes = filtered_nodes # triggers self.update_sel_nodes
    
    @param.depends('node_query', 'filtered_nodes', watch=True)
    def apply_query(self):
        if not self.node_query=='':
            query_nodes = pd.Series(self.node_query.strip().split('\n'))
            str_matches = self.sym_to_index[self.index_col][self.sym_to_index[self.gene_symbol_col].isin(query_nodes[~query_nodes.str.isnumeric()])].unique().tolist()
            int_matches = query_nodes[query_nodes.str.isnumeric()].astype(int).unique().tolist()
            filtered_nodes = self.filtered_nodes[self.filtered_nodes[self.index_col].isin(str_matches+int_matches)]
            
            self.query_found = (filtered_nodes[self.index_col].unique().shape[0], len(np.unique(query_nodes)))
            self.queried_nodes = filtered_nodes # triggers self.update_sel_nodes
        
        else:
            self.query_found = None
            self.queried_nodes = self.filtered_nodes
        
    @param.depends('queried_nodes', watch=True)
    def update_sel_nodes(self):
        sel_nodes = self.annotations.reindex(self.queried_nodes[self.index_col].unique())
        sel_nodes['PPI_SUM_FILT'] = self.compute_PPI_sum(self.queried_nodes)
        self.sel_nodes = sel_nodes
        
    @param.depends('STRINGdb_score', 'sel_nodes', watch=True)
    def update_sel_edges(self):
        in_source = self.edges[self.source_col].isin(self.sel_nodes.index)
        in_target = self.edges[self.target_col].isin(self.sel_nodes.index)
        pass_cutoff = self.edges[self.edge_score_col]>=self.STRINGdb_score
        
        sel_edges = self.edges[in_source & in_target & pass_cutoff].copy()
                  
        PPI_SUM_col = 'PPI_SUM_TOTAL'
        PPI_SUM_filt_col = 'PPI_SUM_FILT'
                  
        PPI_SUM_A = sel_edges[self.source_col].map(self.sel_nodes[PPI_SUM_col])
        PPI_SUM_B = sel_edges[self.target_col].map(self.sel_nodes[PPI_SUM_col])
        sel_edges['min_'+PPI_SUM_col] = np.vstack([PPI_SUM_A, PPI_SUM_B]).min(axis=0)
        
        PPI_SUM_A_filt = sel_edges[self.source_col].map(self.sel_nodes[PPI_SUM_filt_col])
        PPI_SUM_B_filt = sel_edges[self.target_col].map(self.sel_nodes[PPI_SUM_filt_col])
        sel_edges['min_'+PPI_SUM_filt_col] = np.vstack([PPI_SUM_A_filt, PPI_SUM_B_filt]).min(axis=0)
        
        self.sel_edges = sel_edges
        
    @param.depends('max_nodes', 'sel_edges', 'node_display_priority', 'vis_unconnected', watch=True) 
    def update_show_data(self):
        
        node_display_priority = dict(zip(["# PPI observations (all)", "# PPI observations (filtered)"], ['PPI_SUM_TOTAL', 'PPI_SUM_FILT']))[self.node_display_priority]

        if self.vis_unconnected=='Hide':
            temp = []
            for i, x in self.sel_edges.sort_values('min_'+node_display_priority, ascending=False).iterrows():
                temp.extend(x[[self.source_col, self.target_col]].values.tolist())

                if len(np.unique(temp).tolist())>=self.max_nodes:
                    break

            show_nodes = self.sel_nodes[self.sel_nodes.index.isin(np.unique(temp))].reset_index()
        else:
            show_nodes = self.sel_nodes.sort_values(node_display_priority, ascending=False).iloc[:self.max_nodes, :].reset_index()
        
        in_source = self.sel_edges[self.source_col].isin(show_nodes[self.index_col])
        in_target = self.sel_edges[self.target_col].isin(show_nodes[self.index_col])
        show_edges = self.sel_edges[in_source & in_target].copy()
        
        show_nodes['connectivity'] = pd.concat([show_edges.groupby(self.source_col).size(), show_edges.groupby(self.target_col).size()], axis=1).sum(axis=1).reindex(show_nodes[self.index_col]).fillna(0).values
        
        self.param.set_param(show_nodes = show_nodes, show_edges = show_edges) # triggers Network.update_data
