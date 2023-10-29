# MDIC3

## Matrix Decomposition to Infer Cell-Cell Communication

![image1](https://github.com/LYxiaotai/MDIC3/blob/main/Fig1.jpg)

We proposed a new method MDIC3 (Matrix Decomposition to Infer Cell-Cell Communication) to reveal cell-cell communication from cooperation of gene regulatory network (GRN) and matrix decomposition based on scRNA-seq data independent of prior knowledge, and it means the inference of cell-cell communication is no longer limited to special species and does not rely on specific L-R interactions or signaling pathways. MDIC3 can accomplish the inference of intercellular communication networks from individual cells to cell types based on existing cell classification labels. 

## How to use MDIC3

We provide two choices to use MDIC3:
1. You can download the [MDIC3.py](https://github.com/LYxiaotai/MDIC3/tree/main) and [MDIC3_LR.py](https://github.com/LYxiaotai/MDIC3/tree/main) and enter the command line in the Python terminal.
2. You can download the MDIC3 Python package to use the MDIC3 algorithm. MDIC3 Python package can be easily installed:
``` python

pip install MDIC3

```
   
## Inference of cell-cell communication

This tutorial is the example analysis with MDIC3 on a test single-cell gene expression profile that contains 39 genes, 20 cells, 3 cell types.

### 1. Required input data

MDIC3 requires two types of input data:

#### (1). single-cell gene expression data, e.g., [gene_exp.txt](https://github.com/LYxiaotai/MDIC3/tree/main/data/test_data)

* The single-cell gene expression data must be a *.txt file, while each row represents a gene and each column represents a cell. 

* Note that the first row of the single-cell gene expression data must only contains cell names but not contains a header field. An example snapshot of the input single-cell gene expression data format is shown below.

||cell1|cell2|cell3|...|
|-:|:-:|:-:|:-:|:-|
|**Gene1**|3.402|0|8.2|...|
|**Gene2**|1.307|1.012|9.239|...|
|**Gene3**|0|6.051|0|...|
|...|...|...|...|...|

#### (2). scRNAseq metadata, e.g., [cell_label.txt](https://github.com/LYxiaotai/MDIC3/tree/main/data/test_data)

* The scRNAseq metadata must be a *.txt file, while the first column represents cells and the second column represents the corresponding cell labels for the cells in the first column. 

* Note that the first column of the scRNAseq metadata should match exactly with the first row of the single-cell gene expression data. An example snapshot of the input scRNAseq metadata format is shown below.

|cell1|Cell-type1|
|-:|:-|
|**cell2**|**Cell-type1**|
|**cell3**|**Cell-type2**|
|...|...|

We suggest the users to check the two types of input data carefully before running MDIC3. 

We suggest the users to input the single cell gene expression data containing at least 10 cells.

### 2. MDIC3-Inference

The core of MDIC3 is to infer the regulatory relationships among cells based on the regulatory relationships among genes. MDIC3 utilizes GRNs to extract gene regulatory information. We used the [GNIPLR](https://github.com/zyllluck/GNIPLR) algorithm (Gene networks inference based on projection and lagged regression) to infer GRN in our paper. GNIPLR projected gene data twice using the LASSO projection algorithm and the linear projection approximation to produce a linear and monotonous pseudo-time series, and then determined the direction of regulation in combination with lagged regression analyses. You can find more details of GNIPLR in the original article [(doi: 10.1093/bioinformatics/btab099)](https://doi.org/10.1093/bioinformatics/btab099). 
  
#### 2.1 Choice1: Use the "[MDIC3.py](https://github.com/LYxiaotai/MDIC3/tree/main)" to infer cell-cell communications.

It should be noted that the MDIC3 is not limited to GNIPLR. Any tool that infers gene regulatory networks can be used for MDIC3.

Here, we provide two calculation choices, the parameter is shown below:

##### usage1:

*  users can choose to first calculate the GRN using GNIPLR and then infer the cell-cell communications. 

``` python
# Enter the following command line in the Python terminal
python MDIC3.py -exp=scRNA_expression_file -label=cell_label_file -grnchoose='GNIPLR' -process=process_value -step=GRN_calculation_step -out=results_output_fold

# If the users want to save the results of the GRN calculation of GNIPLR, they can do so by adding -grnsave =='TRUE'
python MDIC3.py -exp=scRNA_expression_file -label=cell_label_file -grnchoose='GNIPLR' -process=process_value -step=GRN_calculation_step -grnsave='TRUE' -out=results_output_fold

```

##### usage2:  

* users can also choose to import the GRN calculated by other methods or tools and then infer cell-cell communications using MDIC3.

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
    -grnsave: if -grnsave =='TRUE', the results of the GRN calculated by GNIPLR will be saved to a txt file.
    -grn: if -grnchoose =='other', the user must provide the gene regulatory network adjacency matrix *.txt file.
    -out: the directory to store the MDIC3 results.

* When users choose to calculate the GRN by using GNIPLR, considering that single cell gene expression data always contains a large number of genes and it will take a long time to compute the GRN, we have improved the computational process of GNIPLR through two paths to increasing the computational efficiency of GRN. One path is to use the Python-multiprocessing package process pool, where users should choose the number of processes according to their computer configuration，and the “multiprocessing” package is required. Another path is calculated by blocking the target GRN matrix. When the number of genes is less than 3000, users can set the '-step' parameter to the number of genes or a smaller value than the number of genes. When the number of genes is less than 3000, we suggest users set the '-step' parameter to a value of 3000 or less.

* When users choose to import the GRN results calculated by other methods, users must provide the calculated GRN results, which must be a *.txt file containing only the numerical results of the GRN adjacency matrix, e.g., [GRN.txt](https://github.com/LYxiaotai/MDIC3/tree/main/data/test_data)

* The output of MDIC3 consists of two txt files: 

1) 'celltype_communication.txt' is the result of communication among different cell types, while the first column and the first row represent different cell labels. More specifically, the value in row i and column j indicates the strength of the communication signal sent from the i-th cell type and received by the j-th cell type.

2) 'cellular_communication.txt' is the result of communication among single cells. This is a file that only contains the numeric matrix of the results of communication results among single cells. If the value in row i and column j is greater than 0, it means that the communication signal is sent from the i-th cell and received by the j-th cell. If the value in row i and column j is less than 0, it means that the communication signal is sent from the j-th cell and received by the i-th cell. The absolute value of the value in row i and column j indicates the communication strength between cell i and cell j.

Note that if users choose to calculate the GRN by using GNIPLR, the user can choose whether or not to save the GRN calculation results. If '-grnsave==TRUE', an additional txt file will be generated for saving the GRN adjacency matrix. Note that when the user selects '-grnsave==true', the GRN results file to be saved will be very large, a GRN file containing 20,000 genes may require 10G of storage space. We recommend that users select a storage location with sufficient storage space.


#### There is a simple example below:

The file '[gene_exp.txt](https://github.com/LYxiaotai/MDIC3/tree/main/data/test_data)' is a test single-cell gene expression data that contains 39 genes, 20 cells.

The file '[cell_label.txt](https://github.com/LYxiaotai/MDIC3/tree/main/data/test_data)' is a test scRNAseq metadata that contains 3 cell types corresponding to the 20 cells.

Users can choose two calculation choices to calculate the cell-cell communications results for the 20 cells and the cell type communication results for the 3 cell types. The usages are as follows:

##### usage1: choose to first calculate the GRN using GNIPLR and then infer the cell-cell communications

``` python
# Enter the following command line in the Python terminal
python MDIC3.py -exp=gene_exp.txt -label=cell_label.txt -grnchoose='GNIPLR' -process=2 -step=15 -out=target
# If the users want to save the results of the GRN calculation of GNIPLR, they can do so by adding -grnsave =='TRUE'
python MDIC3.py -exp=gene_exp.txt -label=cell_label.txt -grnchoose='GNIPLR' -process=2 -step=15 -grnsave='TRUE' -out=target

```

The output inluding '[celltype_communication.txt](https://github.com/LYxiaotai/MDIC3/tree/main/data/test_data/result)' and '[cellular_communication.txt](https://github.com/LYxiaotai/MDIC3/tree/main/data/test_data/result)' will be put in your "target" directory. If '-grnsave==TRUE', an additional txt file '[GRN_GNIPLR.txt](https://github.com/LYxiaotai/MDIC3/tree/main/data/test_data/result)' will be generated for saving the GRN adjacency matrix and put in your "target" directory.



##### usage2: choose to import the GRN calculated by other methods or tools and then infer cell-cell communications using MDIC3.

``` python
# The GRN.txt is the GRN adjacency matrix calculation results for the input single-cell gene expression data
# Enter the following command line in the Python terminal
python MDIC3.py -exp=gene_exp.txt -label=cell_label.txt -grnchoose='other' -grn=GRN.txt -out=target
```

The output inluding '[celltype_communication.txt](https://github.com/LYxiaotai/MDIC3/tree/main/data/test_data/result)' and '[cellular_communication.txt](https://github.com/LYxiaotai/MDIC3/tree/main/data/test_data/result)' will be put in your "target" directory. 


#### 2.2 Choice2: Use the MDIC3 Python package to infer cell-cell communications.

* 2.2.1 The MDIC3 Python package also provides two choices to obtain the GRN.
  
  (1) Users can choose to first calculate the GRN using GNIPLR and then infer the cell-cell communication, the function is "lucky.GRN_GNIPLR(AA, gene_exp, step, process)":

      AA：A Python nested list，storing single-cell data expression, with the length of the list equal to the number of genes, the length of each sub-list equal to the number of cells, and the values in the list indicating the expression of each gene on each cell
      gene_exp: a Python dictionary, where each gene name corresponds to a key of the dictionary, and the expression of each gene is stored as a list of numbers in the values of the dictionary.
      step：an integer, the user must select the step size for GNIPLR block calculation.
      process：an integer, the user must select the number of work processes used.

   (2) Users can also choose to import the GRN results calculated by other methods. 

  
* 2.2.2 Users can use the function "lucky.MDIC3_score(AA, GRN, labels, label_index)" to infer the cell-cell communications:

      AA：A Python nested list，storing single-cell data expression, with the length of the list equal to the number of genes, the length of each sub-list equal to the number of cells, and the values in the list indicating the expression of each gene on each cell.
      GRN: A Python numpy array, storing the GRN adjacency matrix.
      labels：A Python list, storing all the cell types, with the length of the list equal to the number of cell types.
      Label_index：a Python dictionary, where each cell type name corresponds to a key of the dictionary, and the index of each cell is stored as a list of numbers in the values of the dictionary.

   The output of this function including the result of communication among single cells and the result of communication among different cell types. 


* 2.2.3 User can save the inferred results using the function "lucky.MDIC3_scoresave(ccc_adjacency, type_adjacency, labels)" (optional)
  
      Parameters：
         ccc_adjacency：the result of communication among single cells 
         type_adjacency: the result of communication among different cell types
         labels：A Python list, storing all the cell types, with the length of the list equal to the number of cell types.


#### There is a simple example below:

##### usage1: choose to first calculate the GRN using GNIPLR and then infer the cell-cell communications

``` python

from MDIC3 import lucky
if __name__ == '__main__':

    AA, gene_exp, cellname = lucky.readexp('gene_exp.txt')
    labels, label_index, label_cell = lucky.readlabel('cell_label.txt')

    # calculate the GRN using GNIPLR
    step = 5        # the user must select the step size for GNIPLR block calculation.
    process = 2     # the user must select the number of work processes used.
    GRN = lucky.GRN_GNIPLR(AA, gene_exp, step, process)

    # Infer the cell-cell communication
    ccc_adjacency, type_adjacency = lucky.MDIC3_score(AA, GRN, labels, label_index)

    # User can save the inferred results using the following command
    lucky.MDIC3_scoresave(ccc_adjacency, type_adjacency, labels)

```

##### usage2: choose to import the GRN calculated by other methods or tools and then infer cell-cell communications using MDIC3.

* When users choose to import the GRN results calculated by other methods. If there is a *.txt file containing only the numerical results of the GRN adjacency matrix, e.g., [GRN.txt](https://github.com/LYxiaotai/MDIC3/tree/main/data/test_data)
  
``` python

from MDIC3 import lucky
if __name__ == '__main__':

    AA, gene_exp, cellname = lucky.readexp('gene_exp.txt')
    labels, label_index, label_cell = lucky.readlabel('cell_label.txt')

    # import the GRN calculated by other methods or tools
    GRN = np.loadtxt(‘GRN.txt')

    # Infer the cell-cell communication
    ccc_adjacency, type_adjacency = lucky.MDIC3_score(AA, GRN, labels, label_index)

    # User can save the inferred results using the following command
    lucky.MDIC3_scoresave(ccc_adjacency, type_adjacency, labels)

```


## Identify L-R pairs from cell-cell communication

The main purpose of applying cell-cell communication analysis is to explain the cell functions through L-R pairs. In our paper, we use a simple method to extract L-R pairs from cell-cell communication.

#### The "[MDIC3_LR.py](https://github.com/LYxiaotai/MDIC3/tree/main)" can be used to extract L-R pairs from cell-cell communication.

#### Required input data

##### 1. single-cell gene expression data, e.g., [gene_exp.txt](https://github.com/LYxiaotai/MDIC3/tree/main/data/test_data)

* The single-cell gene expression data must be a *.txt file, while each row represents a gene and each column represents a cell. 

* Note that the first row of the single-cell gene expression data must only contains cell names but not contains a header field. An example snapshot of the input single-cell gene expression data format is shown below.

||cell1|cell2|cell3|...|
|-:|:-:|:-:|:-:|:-|
|**Gene1**|3.402|0|8.2|...|
|**Gene2**|1.307|1.012|9.239|...|
|**Gene3**|0|6.051|0|...|
|...|...|...|...|...|

##### 2. scRNAseq metadata, e.g., [cell_label.txt](https://github.com/LYxiaotai/MDIC3/tree/main/data/test_data)

* The scRNAseq metadata must be a *.txt file, while the first column represents cells and the second column represents the corresponding cell labels for the cells in the first column. 

* Note that the first column of the scRNAseq metadata should match exactly with the first row of the single-cell gene expression data. An example snapshot of the input scRNAseq metadata format is shown below.

|cell1|Cell-type1|
|-:|:-|
|**cell2**|**Cell-type1**|
|**cell3**|**Cell-type2**|
|...|...|

##### 3. ligand-receptor information, e.g., [LR_human.txt](https://github.com/LYxiaotai/MDIC3/tree/main/data)

* The ligand-receptor information must be a *.txt file, while the first column represents cells and the second column represents the corresponding cell labels for the cells in the first column.

* The file contains only a single column of data with a column name. Each row represents a pair of ligand-receptor pairs, recorded in the format "L - R".

* Note that if a ligand corresponds to a subunit architecture receptor, e.g. IL6 receptors IL6R and IL6ST, it is represented in the format "L - (R1+R2)".

* Considering that [CellChatDB](https://github.com/sqjin/CellChat) L-R database contains both human and mouse ligand-receptor information, we used the human and mouse ligand-receptor genes obtained from CellChatDB to further analyze in the section "Identifying key L-R pairs from cell-cell communication" of our paper. You can download [here](https://github.com/LYxiaotai/MDIC3/tree/main/data).

#### usage:

* Note that the cell type name you write into the command must correspond to the cell type name contained in your scRNAseq metadata file.
* Remember to change the ligand-receptor information according to the type of species you are using.

``` python
# Enter the following command line in the Python terminal
python MDIC3_LR.py -exp='gene_exp.txt' -label='cell_label.txt' -lrdb='LR_human.txt' -ltype='celltype1' -rtype='celltype2' -out=target
```

##### Options and arguments:

    -exp: the input single-cell gene expression data.
    -label: the input scRNAseq metadata.
    -lrdb: the input ligand-receptor information file.
    -ltype: the cell type name that sends cell-cell communication signals during the target communication progress you want to analyze.
    -rtype: the cell type name that receives cell-cell communication signals during the target communication progress you want to analyze.
    -out: the directory to store the L-R identification results. 

* The output of MDIC3_LR.py will be put in your "target" directory. The output file is a txt file named 'target_LR' that includes the L-R pairs involved in the cellular communication process you entered. The output file also contains the Pearson correlation coefficients and significance p-values corresponding to each L-R pair.






## Visualize the MDIC3 inference results

After obtaining the output of MDIC3.py, the 'celltype_communication.txt' can be used to visualize. We provide an R script for plotting the heatmap of the cell type communication results and users can download the script file [here](https://github.com/LYxiaotai/MDIC3/tree/main/Visualize). 

In our paper, we have analyzed the cell type communication of the human lesional skin dataset using MIDC3. We use the inferred results from this dataset as an example to further illustrate how to do visualization for the inference results. Users can download the MDIC3 inferred cell type communications results of the human lesional skin dataset 'celltype_CCC.txt' [here](https://github.com/LYxiaotai/MDIC3/tree/main/Visualize). Users can visualize the results of cell type communications of the human lesional skin dataset by using the following R language code. The heatmap of the cell type communication results will be saved in [plot.pdf](https://github.com/LYxiaotai/MDIC3/tree/main/Visualize). Red-marked parts of the plot indicate significant communication between the two cell types, while the blue-marked parts indicate less significant communication between the two cell types. For example, the square in the first row and the third column of plot.pdf is marked in red, indicating significant communications among inflame_FIB and cDC2 in the predicted results of MDIC3.


```R
library(pheatmap)
library(grid)
df<-read.table("celltype_CCC.txt",header=T,row.names = 1)
df1<-data.matrix(df)
max1<-round(max(df1),2)
min1<-round(min(df1),2)
mean1<-round(mean(df1),2)
step1 = 0.001
bk <- c(seq(min1,mean1,by=step1),seq(mean1+step1,max1,by=step1))
bk1 = seq(min1,mean1,by=step1)
bk2 = seq(mean1+step1,max1,by=step1)
pheatmap(df,
         scale = "none",
         color = c(colorRampPalette(colors = c("Blue","white"))(length(bk1)),colorRampPalette(colors = c("white","red"))(length(bk2))),
         legend = TRUE,
         legend_breaks=seq(min1,max1,mean1-min1),   
         legend_labels=c(min1,mean1,max1),
         breaks=bk,
         cluster_rows=FALSE,
         cluster_cols=FALSE,
         width = 7,height = 6.5,   # users can choose another size of the image to be saved
         filename = 'plot.pdf')
```


