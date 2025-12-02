conda activate scanpy
python
import scanpy as sc
import numpy as np
import os
from tqdm import tqdm
import numpy as np

def preprocess_data(folder_path):
    adata = sc.read_10x_mtx(folder_path, var_names='gene_symbols', cache=False)
    
    sc.external.pp.scrublet(adata, random_state=112)
    
    if 'predicted_doublet' in adata.obs.columns:
        adata = adata[adata.obs['predicted_doublet'] == False].copy()
    else:
        print("Warning: 'predicted_doublet' column not found in adata.obs. Doublet filtering skipped.")
    
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # Adjusted for human mitochondrial gene nomenclature
    
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
    adata = adata[adata.obs.n_genes_by_counts < upper_lim].copy()
    adata = adata[adata.obs.pct_counts_mt < 10].copy()
    
    return adata
data_path = '/mnt/Data/'
cancer_types = os.listdir(data_path)
adata_objects = []

for cancer_type in cancer_types:
    cancer_path = os.path.join(data_path, cancer_type)
    gsm_numbers = os.listdir(cancer_path)
    for gsm_number in gsm_numbers:
        gsm_path = os.path.join(cancer_path, gsm_number)
        sample_types = os.listdir(gsm_path)
        for sample_type in tqdm(sample_types, desc=f"Processing {cancer_type}_{gsm_number}"):
            folder_path = os.path.join(gsm_path, sample_type)
            try:
                adata = preprocess_data(folder_path)
                prefix = f"{cancer_type}_{gsm_number}_{sample_type}_"
                adata.obs_names = [prefix + barcode for barcode in adata.obs_names]
                adata_objects.append(adata)
            except Exception as e:
                print(f"Error processing {folder_path}: {e}")

if adata_objects:
    all_data = sc.concat(adata_objects,join='outer')
    output_path = '/mnt/adata.h5ad'
    all_data.write_h5ad(output_path)
else:
    print("No data processed.")

adata.obs['cell_identifier'] = adata.obs.index
def extract_batch(cell_identifier):
    parts = cell_identifier.split('_')
    return parts[2] if len(parts) > 2 else 'Unknown'
adata.obs['batch'] = adata.obs['cell_identifier'].apply(extract_batch)

batch_counts = adata.obs['batch'].value_counts()

batch_count1 = adata.obs['batch'].nunique()
print(f"Number of unique batches: {batch_count1}")

cell_identifiers = adata.obs['cell_identifier']
prefixes = cell_identifiers.apply(lambda x: x.split('_')[0])
adata.obs['Group'] = prefixes
print(adata.obs[['cell_identifier', 'Group']])
group_counts = adata.obs['Group'].value_counts()
print(group_counts)
data = adata.obs[['batch', 'Group']]
sample_counts = data.groupby('Group')['batch'].nunique()
sample_counts

def extract_cancer_group(cell_identifier):
    parts = cell_identifier.split('_')
    if len(parts) >= 4: 
        cancer_type = parts[0]  
        sample_type = parts[3]  
        return f"{cancer_type}_{sample_type}"
    else:
        return "Unknown" 
adata.obs['Cancer_group'] = adata.obs['cell_identifier'].apply(extract_cancer_group)
print(adata.obs['Cancer_group'].value_counts())
data = adata.obs[['batch', 'Cancer_group']]
sample_counts = data.groupby('Cancer_group')['batch'].nunique()

batches_to_keep = adata.obs['batch'].value_counts()[lambda x: x >= 200].index
adata = adata[adata.obs['batch'].isin(batches_to_keep)].copy()
adata_X_data = adata.X[:10, :10].toarray() if sp.issparse(adata.X) else adata.X[:10, :10]
adata_X_df = pd.DataFrame(adata_X_data, index=adata.obs_names[:10], columns=adata.var_names[:10])
print(adata_X_df)
adata.layers['counts'] = adata.X.copy() 
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
import scipy.sparse
df = pd.DataFrame(adata.X[:5, :5].todense() if scipy.sparse.issparse(adata.X) else adata.X[:5, :5], 
                  index=adata.obs_names[:5], 
                  columns=adata.var_names[:5])
print(df)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata

adata = adata[:, adata.var.highly_variable] 
adata = adata.copy()
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt']) 
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)
plt.savefig(dir+"pca_variance_new.pdf")
adata.write(dir + 'adata_PCA_new.h5ad')
scvi.model.SCVI.setup_anndata(adata, layer = "counts", categorical_covariate_keys=["batch"], continuous_covariate_keys=['pct_counts_mt', 'total_counts'])
model = scvi.model.SCVI(adata)
model.train()
adata.obsm['X_scVI'] = model.get_latent_representation()
adata.layers['scvi_normalized'] = model.get_normalized_expression(library_size = 1e4) 
adata.write(dir + 'scvi_new.h5ad')
sc.pp.neighbors(adata, use_rep = 'X_scVI')
sc.tl.umap(adata)
sc.tl.tsne(adata)
adata.write(dir + 'Before_leiden_data_new.h5ad')
adata = sc.read_h5ad(dir+'Before_leiden_data_new.h5ad')
sc.tl.leiden(adata, resolution = 0.1)
leiden_counts = adata.obs['leiden'].value_counts()
print(leiden_counts) 


colors = ["#E78376","Black", "#81CAD9", "#68B9AA", "#7888AA", "#F0B7A6", "#ABB1C9", "#B1DBD1", 
          "#EBC9DF", "#F292BC", "#b982bc", "#E2E1F0", "#FBB463", "#80B1D3", "#F47F72", 
          "#BDBADB", "#FBF8B4", "#8DD1C6", "#f33b99", "#E7481B", "#1b79af", "#179b73", 
          "#8da0cb", "#bbe173", "#970030", "#3d4a78", "#f6e36d", "#ad5f2c"]
plt.close('all')
plt.rcParams['savefig.dpi'] = 300
plt.figure(figsize=(20, 10))
sc.pl.umap(adata, color='leiden', legend_loc='right margin',#'on data',  
palette=colors, show=False)
plt.tight_layout()
plt.savefig(dir+"umap_leiden2.png",bbox_inches='tight')
adata.write(dir + 'leiden_data_new01.h5ad')
















