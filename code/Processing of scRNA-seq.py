import os
import gzip
import shutil

def gzip_file(file_path):
    """
    Compress a file using gzip format.
    """
    if file_path.endswith('.gz'):
        print(f"Skipping {file_path}, already gzipped.")
        return
    
    try:
        if not os.path.exists(file_path):
            print(f"File {file_path} not found, skipping.")
            return
        

        print(f"File {file_path} exists, proceeding with compression.")
        print(f"File size: {os.path.getsize(file_path)} bytes")
        

        with open(file_path, 'rb') as f_in:
            with gzip.open(file_path + '.gz', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        

        os.remove(file_path)
    
    except FileNotFoundError:

        print(f"File {file_path} was not found during compression, skipping.")
        return
    except Exception as e:

        print(f"Error compressing file {file_path}: {e}")
        return

def compress_files_in_folder(folder_path):
    """
    Compress all .tsv and .mtx files in the folder.
    """
    for root, _, files in os.walk(folder_path):
        for file in files:

            if file.endswith('.gz'):
                continue
            

            if file.endswith('.tsv') or file.endswith('.mtx'):
                file_path = os.path.join(root, file)
                print(f"Compressing {file_path} ...")
                gzip_file(file_path)


base_dir = '/mnt/Data/'


for cancer_type in os.listdir(base_dir):
    cancer_type_path = os.path.join(base_dir, cancer_type)
    
    if os.path.isdir(cancer_type_path):
        for gsm_number in os.listdir(cancer_type_path):
            gsm_path = os.path.join(cancer_type_path, gsm_number)
            
            if os.path.isdir(gsm_path):
                for sample_type in os.listdir(gsm_path):
                    sample_path = os.path.join(gsm_path, sample_type)
                    
                    if os.path.isdir(sample_path):
                        compress_files_in_folder(sample_path)


conda activate scanpy_env
python
import scanpy as sc
import numpy as np
import os
from tqdm import tqdm
import numpy as np  

def preprocess_data(folder_path):
    # Load the data
    adata = sc.read_10x_mtx(folder_path, var_names='gene_symbols', cache=False)
    
    # Predict doublets using scrublet with a fixed random state for reproducibility
    sc.external.pp.scrublet(adata, random_state=112)
    
    # Filter out predicted doublets
    # Ensure 'predicted_doublet' is the correct column name added by scrublet to adata.obs
    if 'predicted_doublet' in adata.obs.columns:
        adata = adata[adata.obs['predicted_doublet'] == False].copy()
    else:
        print("Warning: 'predicted_doublet' column not found in adata.obs. Doublet filtering skipped.")
    
    # Continue with further preprocessing
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    # Identify mitochondrial genes - adjusted for human data where mitochondrial genes start with 'MT-'
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # Adjusted for human mitochondrial gene nomenclature
    
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # Filter based on n_genes_by_counts and pct_counts_mt
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
    output_path = '/mnt/sdb16t/CYJ/Out/adata_qc0827.h5ad'
    all_data.write_h5ad(output_path)
else:
    print("No data processed.")
#######################################################################

python
conda activate scanpy_env
python
import scanpy as sc
import scvi
import os
import math
import itertools
import warnings
import numpy as np
import scipy.sparse as sp
import pandas as pd
import matplotlib.pyplot as plt
adata = sc.read_h5ad('/mnt/adata_qc.h5ad')
dir = '/mnt/Out/'


adata.obs['cell_identifier'] = adata.obs.index
def extract_batch(cell_identifier):
    parts = cell_identifier.split('_')
    return parts[2] if len(parts) > 2 else 'Unknown'
adata.obs['batch'] = adata.obs['cell_identifier'].apply(extract_batch)


batch_counts = adata.obs['batch'].value_counts()

batch_count1 = adata.obs['batch'].nunique()
print(f"Number of unique batches: {batch_count1}")



batches_to_keep = adata.obs['batch'].value_counts()[lambda x: x >= 200].index
adata = adata[adata.obs['batch'].isin(batches_to_keep)].copy()


adata_X_data = adata.X[:10, :10].toarray() if sp.issparse(adata.X) else adata.X[:10, :10]


adata_X_df = pd.DataFrame(adata_X_data, index=adata.obs_names[:10], columns=adata.var_names[:10])



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
plt.savefig(dir+"pca_variance.pdf")
adata.write(dir + 'adata_PCA.h5ad')

adata = sc.read_h5ad('/mnt/adata_PCA.h5ad')
dir = '/mnt/Out/'

####sc.external.pp.bbknn(adata, batch_key= "batch")

scvi.model.SCVI.setup_anndata(adata, layer = "counts", categorical_covariate_keys=["batch"], continuous_covariate_keys=['pct_counts_mt', 'total_counts'])
model = scvi.model.SCVI(adata)
model.train()  #
adata.obsm['X_scVI'] = model.get_latent_representation()
adata.layers['scvi_normalized'] = model.get_normalized_expression(library_size = 1e4)  
adata.write(dir + 'scvi.h5ad')

adata = sc.read_h5ad('/mnt/scvi.h5ad')

sc.pp.neighbors(adata, use_rep = 'X_scVI')
sc.tl.umap(adata)
sc.tl.tsne(adata)
adata.write(dir + 'Before_leiden_data.h5ad')
adata = sc.read_h5ad(dir+'Before_leiden_data.h5ad')
sc.tl.leiden(adata, resolution = 0.1)
#sc.tl.leiden(adata, resolution = 0.2)
leiden_counts = adata.obs['leiden'].value_counts()
print(leiden_counts) 


colors = ["#E78376","Black", "#81CAD9", "#68B9AA", "#7888AA", "#F0B7A6", "#ABB1C9", "#B1DBD1", 
          "#EBC9DF", "#F292BC", "#b982bc", "#E2E1F0", "#FBB463", "#80B1D3", "#F47F72", 
          "#BDBADB", "#FBF8B4", "#8DD1C6", "#f33b99", "#E7481B", "#1b79af", "#179b73", 
          "#8da0cb", "#bbe173", "#970030", "#3d4a78", "#f6e36d", "#ad5f2c"]


dir= '/mnt/Out/'
plt.close('all')

plt.rcParams['savefig.dpi'] = 600
plt.figure(figsize=(20, 10))
sc.pl.tsne(adata, color='leiden', legend_loc='right margin', palette=colors,show=False)
plt.tight_layout()
plt.savefig(dir+"tsne_leiden.png",bbox_inches='tight')

plt.close('all')
plt.rcParams['savefig.dpi'] = 300
plt.figure(figsize=(20, 10))
sc.pl.umap(adata, color='leiden', legend_loc='right margin',#'on data',  
palette=colors, show=False)
plt.tight_layout()
plt.savefig(dir+"umap_leiden2.png",bbox_inches='tight')

adata.write(dir + 'leiden_data.h5ad')
#adata.write(dir + 'leiden_data02.h5ad')
#########################
import matplotlib.pyplot as plt
import scanpy as sc
############################
####################################

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
dir = '/mnt/Out/'
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
all_groups_df = pd.DataFrame()

for group in groups:
    group_data = pd.DataFrame({
        f"{group}_{key}": result[key][group]
        for key in ['names', 'pvals', 'logfoldchanges', 'scores']  
    })
    all_groups_df = pd.concat([all_groups_df, group_data], axis=1)


all_groups_df.to_csv('/mnt/sdb16t/CYJ/Out/all_groups01_data.csv', index=False)



cluster_mapping = {
    "0": "T_NK",
    "1": "Myeloid",
    "2": "Fibroblast",
    "3": "Osteoblast",
    "4": "Epithelium", 
    "5": "B cell",
    "6": "Endothelial",
    "7": "Undefined",
    "8": "Proximal tubule cell",
    "9": "Plasma cell",
    "10": "Epithelium",
    "11": "Epithelium",
    "12": "Lung EpiC",
    "13": "Mast cell",
    "14": "Basal cell",
    "15": "Mix"
}
adata.obs['leiden'] = adata.obs['leiden'].map(cluster_mapping)


colors = (
"#FBB463",
                "#80B1D3",
                "#F47F72",
                "#BDBADB",
                "#D65190",
                "#8DD1C6",
                "#f33b99",
                "#E7481B",
                "#C2D95E",
                "#179b73",
                "#8da0cb",
                'palegreen',
                '#FF7F00',
                "#EBC9DF",
                'dodgerblue',
                "#3d4a78",
                "#1b79af",
                "#f6e36d",
                "#ad5f2c",
                "#FBF8B4"
)
leiden_counts = adata.obs['leiden'].value_counts()
sorted_leiden = leiden_counts.index.tolist()
adata.obs['leiden'] = pd.Categorical(adata.obs['leiden'], categories=sorted_leiden, ordered=True)
################
cell_types = adata.obs['leiden'].unique()

palette = dict(zip(cell_types, colors))
plt.rcParams['savefig.dpi'] = 600
fig = sc.pl.umap(
    adata,
    color='leiden', 
    palette=colors,  
    legend_loc='right margin',  
    title='UMAP',
    show=False,  
    return_fig=True  
)
fig.set_size_inches(5, 5)
plt.savefig(dir + "umap.png", bbox_inches='tight')
plt.close()

#####################################
#####################################
#####################################


##############################

marker_genes = [
 "CD3D","KLRD1","GZMA",
"CD14","MS4A7","APOE",
"KRT18","KRT19","EPCAM",
"COL1A1","DCN","POSTN",
"IBSP","CLEC11A","ALPL",
"CD79A","CD79B","MS4A1",
"PLVAP","VWF","CDH5",
"SLC3A1","SLC22A6","SLC17A4",
"IGKC","MZB1","JCHAIN",
"SFTPC", "SFTPB", "NAPSA", 
"TPSAB1","TPSB2","CPA3",
"KRT14","KRT5","TP63"
]
plt.rcParams.update({'font.size': 12})
dotplot = sc.pl.dotplot(
  adata, 
  marker_genes, 
 cmap="plasma", 
  groupby='leiden', 
)  
fig = plt.gcf()
fig.set_size_inches(13, 5)  
dir="/mnt/Out/"
plt.savefig(dir+"dotplot.pdf", bbox_inches='tight')


########################################
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap



custom_colors = ["#fad1e0", "#E59CC4", "#984EA3"]
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", custom_colors)


plt.rcParams.update({'font.size': 12})

dotplot = sc.pl.dotplot(
    adata, 
    marker_genes, 
    cmap=custom_cmap,  
    groupby='leiden', 
    show=False  
)

fig = plt.gcf()
fig.set_size_inches(13, 5)


plt.savefig(dir + "dotplot1.pdf", bbox_inches='tight')
plt.close()
