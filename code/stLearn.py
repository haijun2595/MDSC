
conda install -c conda-forge stlearn
############################################################
conda activate stlearn 
python
import stlearn as st
import scanpy as sc
import numpy as np
import os
import random
import pandas as pd
import matplotlib.pyplot as plt
random.seed(1000000)
np.random.seed(1000000)
data_dir = "/mnt/GSM7974888/"
data = st.Read10X(data_dir) 
data.var_names_make_unique() 
print(data.uns["spatial"].keys()) 
#dict_keys(['PT_ST_Visium'])
st.add.image(adata=data,
             imgpath=data_dir+"spatial/tissue_hires_image.png",
             library_id="Sample_5687-AJ-S2-D_GTCGCGACAA-TGGCAGATTG", visium=True)


#st.pp.filter_genes(data, min_cells=3)   
#sc.pp.filter_genes(data, min_counts=10)

data.obs['nFeature_Spatial'] = np.sum(data.X > 0, axis=1)  
data.obs['nCount_Spatial'] = np.sum(data.X, axis=1) 
data.obs['percent_mt'] = np.sum(data[:, data.var_names.str.startswith('MT-')].X, axis=1).A1 / np.sum(data.X, axis=1).A1 * 100
data.obs['percent_rb'] = np.sum(data[:, data.var_names.str.startswith('RP[LS]')].X, axis=1).A1 / np.sum(data.X, axis=1).A1 * 100
data = data[(data.obs['nFeature_Spatial'] >= minFeature) & 
                         (data.obs['percent_mt'] <= maxMT) &
                         (data.obs['nCount_Spatial'] >= minCount), :]


st.pp.normalize_total(data)
 



lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='human')
print(len(lrs))

np.random.seed(1000000) 

st.tl.cci.run(data,
                        lrs, min_spots = 20, 
                        distance=None, 
                        n_pairs=10000, 
                         n_cpus=80)



lr_info = data.uns['lr_summary']
print('\n', lr_info) 
#1192

st.tl.cci.adj_pvals(data, correct_axis='spot', pval_adj_cutoff=0.05, adj_method='fdr_bh')


# Showing the rankings of the LR from a global and local perspective.# Ranking based on number of significant hotspots.
st.pl.lr_summary(data, n_top=500)

plt.savefig(data_dir+'s2.lr_summary500.pdf',dpi=600)
plt.close()
st.pl.lr_summary(data, n_top=50,figsize=(10,3))


st.pl.lr_summary(data, n_top=30, figsize=(9,6))
plt.savefig(data_dir+'s2.lr_summary30.pdf',dpi=600)
plt.close()

st.pl.lr_diagnostics(data, figsize=(6,2.5))
plt.savefig(data_dir+'s2.lr_diagnostics.pdf',dpi=600)
plt.close()

st.pl.lr_n_spots(data, n_top=30, figsize=(20,12),max_text=100)
plt.savefig(data_dir+'s2.lr_n_spots30.pdf',dpi=600)
plt.close()

st.pl.lr_n_spots(data, n_top=500, figsize=(11, 4), max_text=100)
plt.savefig(data_dir+'s2.lr_n_spots500.pdf',dpi=600)
plt.close()


r_path = "/home/ps/R-4.0.5/doc/RESOURCES"
st.tl.cci.run_lr_go(data, r_path)
st.pl.lr_go(data, lr_text_fp={'weight': 'bold', 'size': 10}, rot=15,figsize=(12,3.65), n_top=15, show=False)

best_lr = data.uns['lr_summary'].index.values[0]


stats = ['lr_scores', 'p_vals', 'p_adjs', '-log10(p_adjs)']
fig, axes = plt.subplots(ncols=len(stats), figsize=(16,6))
for i, stat in enumerate(stats):    
           st.pl.lr_result_plot(data, use_result=stat,size=5, use_lr=best_lr, show_color_bar=True, ax=axes[i])
           axes[i].set_title(f'{best_lr} {stat}')

plt.savefig(data_dir+'s4.lr_result_.pdf',dpi=600)
plt.close()

fig, axes = plt.subplots(ncols=2, figsize=(8,6))  

st.pl.lr_result_plot(data, 
         use_result='-log10(p_adjs)',
          use_lr=best_lr, 
           show_color_bar=False,
           ax=axes[0])


st.pl.lr_result_plot(data, 
               use_result='lr_sig_scores',
                use_lr=best_lr, 
                  size=2,
                   show_color_bar=False,
                          ax=axes[1])
axes[0].set_title(f'{best_lr} -log10(p_adjs)')
axes[1].set_title(f'{best_lr} lr_sig_scores')
plt.savefig(data_dir+'s4.lr_result_plot.pdf',dpi=600)
plt.close()

import pandas as pd


lr_summary_index = data.uns['lr_summary'].index.values
df = pd.DataFrame(lr_summary_index, columns=['Interaction'])
df.to_csv(data_dir+'lr_summary_interactions.csv', index=False)




selected_lrs = ["COL1A1_CD44","FN1_CD44","COL1A2_CD44"]
selected_lrs = ["FN1_CD44"]
stats = ['lr_scores', 'p_adjs', '-log10(p_adjs)','lr_sig_scores'] #'p_vals',
for lr in selected_lrs:
    fig, axes = plt.subplots(ncols=len(stats), figsize=(4 * len(stats), 6))  
    for j, stat in enumerate(stats):
           st.pl.lr_result_plot(data, use_result=stat,use_lr=lr,size=9,show_color_bar=True,ax=axes[j])
           axes[j].set_title(f'{lr} {stat}')
           handles, labels = axes[j].get_legend_handles_labels()        
           if handles:            
                      axes[j].legend(handles=handles, labels=labels, loc='best')
plt.tight_layout() 
plt.savefig(f'{data_dir}{lr}_result_plots.pdf')
plt.close()  

