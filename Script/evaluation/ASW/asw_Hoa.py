def silhouette_coeff_ASW(adata, method_use='raw',save_dir='', save_fn='', percent_extract=0.8):
    random.seed(0)
    asw_fscore = []
    asw_bn = []
    asw_bn_sub = []
    asw_ctn = [] 
    iters = []
    for i in range(20):
        iters.append('iteration_'+str(i+1))
        rand_cidx = np.random.choice(adata.obs_names, size=int(len(adata.obs_names) * percent_extract), replace=False)
        adata_ext = adata[rand_cidx,:]
        asw_batch = silhouette_score(adata_ext.X, adata_ext.obs['batch'])
        asw_celltype = silhouette_score(adata_ext.X, adata_ext.obs['cell_type'])
        min_val = -1
        max_val = 1
        asw_batch_norm = (asw_batch - min_val) / (max_val - min_val)
        asw_celltype_norm = (asw_celltype - min_val) / (max_val - min_val)
        
        fscoreASW = (2 * (1 - asw_batch_norm)*(asw_celltype_norm))/(1 - asw_batch_norm + asw_celltype_norm)
        asw_fscore.append(fscoreASW)
        asw_bn.append(asw_batch_norm)
        asw_bn_sub.append(1-asw_batch_norm)
        asw_ctn.append(asw_celltype_norm)
    

    df = pd.DataFrame({'asw_batch_norm':asw_bn, 'asw_batch_norm_sub': asw_bn_sub,
                       'asw_celltype_norm': asw_ctn, 'fscore':asw_fscore,
                       'method_use':np.repeat(method_use, len(asw_fscore))})
    df.to_csv(save_dir + save_fn + '.csv')
    print('Save output of pca in: ',save_dir)
    return df
        