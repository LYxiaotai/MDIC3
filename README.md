# MDIC3

## Matrix Decomposition to Infer Cell-Cell Communication

![image](https://github.com/LYxiaotai/MDIC3/blob/main/Fig1_20221017.tif)

MDIC3 (Matrix Decomposition to Infer Cell-Cell Communication) is a method for inferring cell-cell communication.

The "MDIC3_CCC.py" can be used to construct the cell-cell communication network.

The "MDIC3_type.py" can be used to construct the communication network among cell types.

There is a simple example below:

The file 'gene_exp.txt' is a single-cell gene expression profile that contains 39 genes, 20 cells.

The file 'GRN.txt' is gene regulatory network adjacency matrix among the 39 genes in 'gene_exp.txt'.

The file 'cell_label.txt' contains cell types informations for the 20 cells in 'gene_exp.txt'.

Step1: Input the file 'gene_exp.txt' and 'GRN.txt';

Step2: Run 'MDIC3_CCC.py';

Step3: Output the result file 'CCC_results.txt'.

The output file 'CCC_results.txt' is cell-cell communication network adjacency matrix among the 20 cells in 'gene_exp.txt'.

Step4: Input the file 'cell_label.txt'and 'CCC_results.txt';

Step5: Run 'MDIC3_type.py';

Step6: Output the result file 'type_results.txt'.

The output file 'type_results.txt' is communication network among different cell types in 'cell_label.txt' and 'gene_exp.txt'.

