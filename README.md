# MDIC3

## Matrix Decomposition to Infer Cell-Cell Communication

![image1](https://github.com/LYxiaotai/MDIC3/blob/main/Fig1_20221017.jpg)

We proposed a new method MDIC3 (Matrix Decomposition to Infer Cell-Cell Communication) to reveal cell-cell communication from cooperation of gene regulatory network and matrix decomposition based on scRNA-seq data independent of prior knowledge, and it means the inference of cell-cell communication is no longer limited to special species and does not rely on specific L-R interactions or signaling pathways. MDIC3 can accomplish the inference of intercellular communication networks from individual cells to cell types based on existing cell classification labels. 

## How to use MDIC3

This tutorial is the example analysis with MDIC3 on a test single-cell gene expression profile that contains 39 genes, 20 cells, 3 cell types.

### Required input data

MDIC3 requires two types of input data:

#### 1. single-cell gene expression data, e.g., gene_exp.txt in ....

  ● The single-cell gene expression data must be a *.txt file, while each row represents a gene and each column represents a cell. 

  ● Note that the first row of the single-cell gene expression data must only contains cell names but not contains a header field. 

#### 2. scRNAseq metadata, e.g., cell_label.txt in ...

  ● The scRNAseq metadata must be a *.txt file, while the first column represents cells and the second column represents the corresponding cell labels for the cells in the first column. 

  ● Note that the first column of the scRNAseq metadata should match exactly with the first row of the single-cell gene expression data.

We suggest the users to check the two types of input data carefully before running MDIC3. 

We suggest the users to input the single cell gene expression data containing at least 10 cells.

### MDIC3-Inference

The core of MDIC3 is to infer the regulatory relationships among cells based on the regulatory relationships among genes. MDIC3 utilizes GRNs to extract gene regulatory information. We used the GNIPLR algorithm (Gene networks inference based on projection and lagged regression) to infer GRN in our paper. GNIPLR projected gene data twice using the LASSO projection algorithm and the linear projection approximation to produce a linear and monotonous pseudo-time series, and then determined the direction of regulation in combination with lagged regression analyses. You can find more details of GNIPLR in the original article (doi: 10.1093/bioinformatics/btab099). It should be noted that the MDIC3 is not limited to GNIPLR. Any tool that infers gene regulatory networks can be used for MDIC3.

The "MDIC3.py" can be used to infer cell-cell communications, and the parameter is shown below:

usage1:

    python MDIC3.py -exp=scRNA_expression_file -label=cell_label_file -grnchoose='other' -grn=grn_file -out=results_output_fold

usage2:

    python MDIC3.py -exp=scRNA_expression_file -label=cell_label_file -grnchoose='GNIPLR' -process= process_value -step=GRN_calculation_step -out=results_output_fold
    
Options and arguments:

    -exp: the input single-cell expression profile
    -label: cell labels for each cell
    -grnchosse: availability of gene regulatory networks
    -grn: if grnchoose =='other', the user must provide the gene regulatory network adjacency matrix txt file.
    -process: if grnchoose =='GNIPLR', the user must select the number of work processes used.
    -step: if grnchoose =='GNIPLR', the user must select the step size for GNIPLR block calculation.
    -out: the directory to store the MDIC3 results    

There is a simple example below:

The file 'gene_exp.txt' is a test single-cell gene expression profile that contains 39 genes, 20 cells.

The file 'cell_label.txt' contains cell types informations for the 20 cells in 'gene_exp.txt'.

Step1: Input the file 'gene_exp.txt' and 'GRN.txt';

Step2: Run 'MDIC3_CCC.py';

Step3: Output the result file 'CCC_results.txt'.

The output file 'CCC_results.txt' is cell-cell communication network adjacency matrix among the 20 cells in 'gene_exp.txt'.

Step4: Input the file 'cell_label.txt'and 'CCC_results.txt';

Step5: Run 'MDIC3_type.py';

Step6: Output the result file 'type_results.txt'.

The output file 'type_results.txt' is communication network among different cell types in 'cell_label.txt' and 'gene_exp.txt'.

