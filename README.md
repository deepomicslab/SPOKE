# SPOKE
A tool for annotating how cells spatially SPOKE to each other.

---
## Envrionment

<!--  Please see requirements.txt for Python environment (Python=3.10.9).--> 
<!--  For R script, ggplot2, igraph, dplyr, and tidyr are required.--> 
For python environment, we offer packed environment file spoke_env.tar.gz at https://yunpan.tongji.edu.cn/link/ARB8AB3BAF2A5A461FA6DB533C3EDC92DE 
You can move it to the anaconda path (like <home>/anaconda/envs/) and unzip the file using:
```
tar -xzf SPOKE21_env.tar.gz
source activate
```

You can also install the packages by yourself:
```
conda create -n spoke python=3.9.10
conda activate spoke
pip install squidpy==1.3.1
pip install pulp==2.7.0
```

## Pipeline

![Pipeline of SPOKE](pipeline_graph7.png)

The processing procedures SPOKE contains are:
- Construction of interaction interfaces (to capture interactions between cell pairs) for cells.
- Construction of the cell graph and the interface graph.
- Construction of motifs (to capture interactions among three cells).
- Finding edging cells (cells participating distinct interface clusters simultaneously).

The analyses SPOKE contains are:
- Construction of the abstract graph representation.
- Enrichment of LR pairs or pathways for interface clusters
- Enrichment of genes or pathways for edging cells.
- Enrichment of LR pairs or pathways for motifs.

## Quick start
For example, we provide the human melanoma data in the example_data folder. 
We will show how to finish the process the input data and finish the analyses.

### Preparation
With the input anndata file, you need to finish the several steps in the following order.
We provide the melanoma data (mentioned in paper) in the example_data folder (named as mela_cell_adata.h5ad).
In the following steps, we take the melanoma data as example.

#### 1. Construction of interface profile
1.a if using h5ad format input:
Note that if there are too many all-zero rows, the all-zeros rows should be deleted in advance.
Given the h5ad file, you can run command like:
```
cd SPOKE/scripts
python "build_infc.py" --mode adata --adata "../example_data/mela_cell_adata.h5ad" --x_col x --y_col y --cell_id_col barcode --gene_id_col gene_id --gene_symbol_col gene_name
```

1.b if using csv format input

Given csv files of cell-gene count matrix, cell meta, and gene meta (optional),
you can run command like:
```
cd SPOKE/scripts
python "build_infc.py" --mode csv --count "../example_data/count.csv" --cellmeta "../example_data/cell_meta.csv" --genemeta "../example_data/gene_meta.csv" --x_col x --y_col y --cell_id_col barcode --gene_id_col gene_id --gene_symbol_col gene_name
``` 
Although we provide such optional way in this step, it is necessary to generate related
.h5ad file for step 2.2.

The constructed interface data are saved under '../tmp' path.
The h5ad file of interface data for the melanoma case is provided in example_data/example folder (named as idata.h5ad).
#### 2. Iterative construction of cell graph and interface graph
2.1 Generate embeddings of the interface profile using STAGUE
We recommend you create additional anaconda environment for running STAGUE to avoid environmental conflict issues (To get further details can be found at https://github.com/deepomicslab/STAGUE.).
You can use command like:
```
cd STAGUE
python main.py --adata_file <home>/SPOKE_/tmp/idata.h5ad --output_dir ./result/V1 --n_clusters 7 --adj_weight 0.3 --k 25 --clu_model kmeans
```
Note that if there are too many interfaces with all-zero signaling value, you need to exclude these rows before generating the embeddings. Too many all-zero rows may badly affect the performance of further analysis.


The LR database file is in the 'data' folder. 

2.2 Iteratively generate cell graph and interface graph. The .h5ad files of cell and interface 
are needed as input (after previous steps, you may have obtained the .h5ad files). You can use command like:
```
cd SPOKE/scripts
python iter_graphs.py --cell_adata ../example_data/mela_cell_adata.h5ad --infc_adata ../example_data/example/idata.h5ad
```
Parameters:
```
--cell_adata: Type: str; Required. Specifies the file path of the cell data in AnnData format.
--infc_adata: Type: str; Required. Specifies the file path of the interface data in AnnData format.
--thres: Type: float; Optional; Default: 0.005. Sets the threshold for the iteration process.
--max_iter: Type: int; Optional; Default: 5. Specifies the maximum number of iterations allowed.
--gamma_C: Type: float; Optional; Default: 0.95. Represents the proportion parameter for the cell graph.
--gamma_I: Type: float; Optional; Default: 0.05. Represents the proportion parameter for the interface graph.
--K: Type: int; Optional; Default: 3. Defines the number of nearest neighbors to consider.
--output: Type: str; Optional; Default: "..". Specifies the output directory where the results will be saved.
```
The iteration process may take some time. Please wait patiently.
The cell graph and interface graph are saved as adjacency matrices in .txt format 
(named as 'S_C.txt' and 'S_I.txt' in the 'example_data/example' folder).

#### 3. Abstract cell or interface clusters
For abstracting we use the clustering algorithm TEC-U. 
Further details can be found at https://github.com/keyqing5/TEC.
Note that it is okay to use any other clustering algorithm, such as Louvain and Leiden.
But please pay attention to the format of output for further computation. 
The cluster labels should be stored in the h5ad file. 
Specifically, the clustering result for cell data cell.h5ad file (read as variable cell_adata in Python) 
should be stored in the 'optimal_cluster' column of cell_adata.obs (i.e., cell_adata.obs['optimal_cluster'])

We save the clustering result of interface data for melanoma case in the example_data/example folder.
To integrate the abstracting result, you can run command like:
```
cd SPOKE/scripts
python 
```

### Analyses

#### Construction of the abstract graph representation
Before running analyses, gseapy package need to be installed by:
```pip install gseapy```.
You can get the adjacent matrix of the abstract graph with the following command:



#### Analysis of LR pairs or pathways for interface clusters and edging cells
We provide the script to run analyses for interface clusters and edging cells together.
You can run command like:
```
cd SPOKE/scripts
python analyz_interface_n_boundary_cell_v2.py <idata name> <adata name> <cluster result file> <save root dir> <cell type key(optional)> <save gene exp map flag(optional)>
```
For the melanoma data, you can run command like:
```
python analyz_interface_n_boundary_cell_v2.py ../example_data/example/idata.h5ad ../example_data/mela_cell_adata.h5ad ../example_data/example/infc_abstract_result.csv "../example_data/example/infc_bc/" cell_type
```
Warning messages of scanpy package can be ignored.

#### Enrichment of LR pairs or pathways for motifs.
Before computation, you need to clone the Network-Motif under SPOKE folder.
```
cd SPOKE
git clone https://github.com/gtremper/Network-Motif.git
```

You can run command like:
```
cd SPOKE/scripts
python analyz_color_motif.py <adata name> <idata name> <cell type key> <save dump files dir> <delete adj file or not[optional]>
```
For the melanoma data, you can run command like:
```
python analyz_color_motif.py ../example_data/mela_cell_adata.h5ad ../example_data/example/idata.h5ad cell_type /root/example_data/example/cm_dump/
```
Note that the 'save dump files dir' parameter need to be an absolute path. This is due to the settings of Network-Motif.
The computation of Network-Motif may need some time. Please wait patiently.

## Input
There are two choices for input:
1. The h5ad file containing cell expression matrix and spatial coordinates.
2. The csv files of cell gene expression count and cell meta information (with spatial coordinates). Gene meta file is optional.

Example is provided in the example_data folder.

## Contact

Feel free to open an issue in Github or contact **liangqs@tongji.edu.cn** and **lingxi.chen@cityu.edu.hk** if you have any problem in using SPOKE.
