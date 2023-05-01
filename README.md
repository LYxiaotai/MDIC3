# MDIC3

## Matrix Decomposition to Infer Cell-Cell Communication

![image1](https://github.com/LYxiaotai/MDIC3/blob/main/Fig1_20221017.jpg)

We proposed a new method MDIC3 (Matrix Decomposition to Infer Cell-Cell Communication) to reveal cell-cell communication from cooperation of gene regulatory network (GRN) and matrix decomposition based on scRNA-seq data independent of prior knowledge, and it means the inference of cell-cell communication is no longer limited to special species and does not rely on specific L-R interactions or signaling pathways. MDIC3 can accomplish the inference of intercellular communication networks from individual cells to cell types based on existing cell classification labels. 

## How to use MDIC3

This tutorial is the example analysis with MDIC3 on a test single-cell gene expression profile that contains 39 genes, 20 cells, 3 cell types.

### Required input data

MDIC3 requires two types of input data:

#### 1. single-cell gene expression data, e.g., gene_exp.txt in ....

* The single-cell gene expression data must be a *.txt file, while each row represents a gene and each column represents a cell. 

* Note that the first row of the single-cell gene expression data must only contains cell names but not contains a header field. 

#### 2. scRNAseq metadata, e.g., cell_label.txt in ...

* The scRNAseq metadata must be a *.txt file, while the first column represents cells and the second column represents the corresponding cell labels for the cells in the first column. 

* Note that the first column of the scRNAseq metadata should match exactly with the first row of the single-cell gene expression data.

We suggest the users to check the two types of input data carefully before running MDIC3. 

We suggest the users to input the single cell gene expression data containing at least 10 cells.

### MDIC3-Inference

The core of MDIC3 is to infer the regulatory relationships among cells based on the regulatory relationships among genes. MDIC3 utilizes GRNs to extract gene regulatory information. We used the [GNIPLR](https://github.com/zyllluck/GNIPLR) algorithm (Gene networks inference based on projection and lagged regression) to infer GRN in our paper. GNIPLR projected gene data twice using the LASSO projection algorithm and the linear projection approximation to produce a linear and monotonous pseudo-time series, and then determined the direction of regulation in combination with lagged regression analyses. You can find more details of GNIPLR in the original article [(doi: 10.1093/bioinformatics/btab099)](https://doi.org/10.1093/bioinformatics/btab099). 
  
#### The "MDIC3.py" can be used to infer cell-cell communications.

It should be noted that the MDIC3 is not limited to GNIPLR. Any tool that infers gene regulatory networks can be used for MDIC3.

Here, we provide two calculation pathways, the parameter is shown below:

##### usage1:

* user can choose to first calculate the GRN using GNIPLR and then calculate the cell-cell communication results. 

``` python
# Enter the following command line in the Python terminal
python MDIC3.py -exp=scRNA_expression_file -label=cell_label_file -grnchoose='GNIPLR' -process= process_value -step=GRN_calculation_step -out=results_output_fold
```

##### usage2:  

* user can also choose to import the GRN results calculated by other methods to calculate cell-cell communication results.

``` python
# Enter the following command line in the Python terminal
python MDIC3.py -exp=scRNA_expression_file -label=cell_label_file -grnchoose='other' -grn=grn_file -out=results_output_fold
```

    
##### Options and arguments:

    -exp: the input single-cell gene expression data.
    -label: the input scRNAseq metadata.
    -grnchoose: availability of gene regulatory networks.
    -process: if -grnchoose =='GNIPLR', the user must select the number of work processes used.
    -step: if -grnchoose =='GNIPLR', the user must select the step size for GNIPLR block calculation.
    -grn: if -grnchoose =='other', the user must provide the gene regulatory network adjacency matrix *.txt file.
    -out: the directory to store the MDIC3 results.

* When users choose to calculate the GRN by using GNIPLR, considering that single cell gene expression data always contains a large number of genes and it will take a long time to compute the GRN, we have improved the computational process of GNIPLR through two paths to increasing the computational efficiency of GRN. One path is to use the Python-multiprocessing package process pool, where users should choose the number of processes according to their computer configuration， and the “multiprocessing” package is required. Another path is calculated by blocking the target GRN matrix. When the number of genes is less than 3000, users can set the '-step' parameter to the number of genes or a smaller value than the number of genes. When the number of genes is less than 3000, we suggest users set the '-step' parameter to a value of 3000 or less.

* When users choose to import the GRN results calculated by other methods, users must provide the calculated GRN results, which must be a *.txt file containing only the numerical results of the GRN adjacency matrix, e.g., GRN.txt in ...

#### There is a simple example below:

The file 'gene_exp.txt' is a test single-cell gene expression data that contains 39 genes, 20 cells.

The file 'cell_label.txt' is a test scRNAseq metadata that contains 3 cell types corresponding to the 20 cells.

We use two calculation pathways to calculate the cell-cell communications results for the 20 cells and the cell type communication results for the 3 cell types, and the results will be put in "target" directory:

##### usage1: 

* choose to first calculate the GRN using GNIPLR and then calculate the cell-cell communication results.


``` python
# Enter the following command line in the Python terminal
python MDIC3.py -exp=gene_exp.txt -label=cell_label.txt -grnchoose='GNIPLR' -process=2 -step=39 -out=target

```

##### usage2: 

* choose to import the GRN results calculated by other methods to calculate cell-cell communication results.


``` python
# The GRN.txt is the GRN adjacency matrix calculation results for the input single-cell gene expression data
# Enter the following command line in the Python terminal
python MDIC3.py -exp=gene_exp.txt -label=cell_label.txt -grnchoose='other' -grn=GRN.txt -out=target
```


## Visualize the MDIC3 inference results

