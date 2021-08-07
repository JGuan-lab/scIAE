# scIAE: an integrative autoencoder-based ensemble classification framework for single-cell RNA-seq data </br> 
## 1. Introduction  
  scIAE is an integrative autoencoder-based ensemble classification framework for single-cell RNA-seq data to identify cell type and predict disease status.
  
  Giving gene expression matrices and labels (cell type annotation or disease status) of training set and testing set, 
 
  scIAE correspond to the following paper:
  
  Yin, Q., Wang, Y., Guan, J., Ji, G. . scIAE: an integrative autoencoder-based ensemble classification framework for single-cell RNA-seq data.
  
## 2. Installation
Depends: 

       R (>= 4.0.4)    
Requiremens: 

      library('edgeR')     
      library('igraph')     
Get CtsDGM from github:  

      install.packages('devtools')    
      devtools::install_github("JGuan-lab/CtsDGM")    
## 3. Quick start
Import the package:

      library('CtsDGM')  
Prepare data: 

The datasets analyzed in the paper are available at: https://doi.org/10.5281/zenodo.4288928
  
      load("exp_count.rda")  #a matrix or data frame of read counts with row and column names, rows denoting genes and column denoting cells.
      load("celltype_annotation.rda")  #a list of cell type annotation of cells. 
      load("ref_network.rda") #a referenced gene-gene interaction network, a two-column matrix, character or numeric, each row defining one edge.
      load("ASD_genes.rda") #a list of disease-associated genes.
        
       
Run the functions:

      gene_specificity = calc_specificity(datExpr = exp_count,cell_type = cell_type)  
      gene_score = calc_score(specificity = gene_specificity,cell_type = cell_type)  
      net_specific = constr_specific_net(net = refnet,genescore = gene_score) 
      Cts_DGM = constr_disease_module(Cts_list = net_specific,disease_gene,padjCutoff = 0.1,returnAllmodules = FALSE,itera = 1000)
## 4. Step by step
### 4.1 Library and load data
      library('CtsDGM')  

The datasets analyzed in the paper are available at: https://doi.org/10.5281/zenodo.4288928
      
      load("exp_count.rda")  
      load("celltype_annotation.rda")   
      load("ref_network.rda")
      load("ASD_genes.rda")

      > exp_count[1:5,1:5]  
            F1S4_160106_001_C01 F1S4_160106_001_G01 F1S4_160106_001_H01 F1S4_160106_002_F01 F1S4_160106_002_H01  
     1                       0                   0                   0                   0                   1  
     29974                   0                   0                   0                   0                   0  
     2                       0                 185                  54                  93                 256  
     53947                   0                   0                   1                   0                   0  
     51146                   0                   0                   0                   0                   0  
      > head(cell_type)  
      [1] "Glutamatergic neuron" "Glutamatergic neuron" "Glutamatergic neuron" "Glutamatergic neuron"  
      [5] "Glutamatergic neuron" "Glutamatergic neuron"  
      > head(disease_gene)  
      [1] 10349 10347  1636    43    60 51412  
      > head(refnet)  
           V1     V2
      1243  1 253982  
      2241  1   4782  
      3194  1  80726  
      4030 10  10290  
      4065 10  10788  
      4089 10  1104  
### 4.2 Calculate cell type specificity of genes
`calc_specificity()` calculates specificity for each cell type and each gene. Its inputs include a expression dataset with rows denoting genes and columns denoting cells, and a cell type annotation infomation for cells.

       gene_specificity = calc_specificity(datExpr = exp_count,cell_type = cell_type) 
       
       > head(gene_specificity)  
             Glutamatergic neuron GABAergic interneuron    astrocyte Oligodendrocyte precursor cell Oligodendrocyte  
       1              0.153923981           0.092971540 0.0038447204                   0.0004915768      6.49671348  
       29974          1.162882650           0.650760445 0.4236547673                   0.2294388593      0.24206465  
       2              0.006847725           0.006332535 0.0155710509                   0.0199735027      0.02583191  
       53947          0.048084097           0.057328215 0.1946284907                   0.0414344054      0.48405174  
       51146          0.006450705           0.001296425 0.0004588732                   0.0006520376      0.01061721  
       8086           0.294264024           0.293470977 0.8562424608                   0.4765065715      1.16789349  
                     No  Microglial endothelial  
       1     0.12440888         NaN         NaN  
       29974 0.85993200  0.07287314   0.0000000  
       2     0.01154272  5.13526653   0.1947319  
       53947 0.15760180  0.15616393   2.0658948  
       51146 0.01649277 60.63264702   0.0000000  
       8086  0.32362597  0.82168237   0.4923334  
### 4.3 Calculate cell type score of genes 
`calc_score()` calculates cell type score for each gene by using the cell type specificity obtained from `calc_specificity()`.

       gene_score = calc_score(specificity = gene_specificity,cell_type = cell_type)
       
       > head(gene_score)  
            Glutamatergic neuron GABAergic interneuron   astrocyte Oligodendrocyte precursor cell Oligodendrocyte  
      29974            1.6187487             0.6199845  0.17707267                    -0.20169607      -0.1770727  
      2               -0.1893733            -0.1983039 -0.03815748                     0.03815748       0.1397109  
      53947           -0.5132814            -0.4696703  0.17807302                    -0.54465269       1.5434888  
      51146            0.2244436            -0.2244436 -0.29738620                    -0.28056347       0.5873050  
      8086            -0.3699266            -0.3714694  0.72333810                    -0.01539462       1.3296195  
      65985            1.0437389             0.3832839 -0.31043785                    -0.22533974       0.2253397  
                     No    Microglial endothelial  
      29974  1.027920472 -5.070377e-01 -0.64915822  
      2     -0.107987208  8.870997e+01  3.06753250  
      53947  0.003391737 -3.391737e-03  9.00617034  
      51146  1.099008531  5.280172e+03 -0.33734955  
      8086  -0.312806345  6.561054e-01  0.01539462  
      65985  1.925996367 -8.750958e-01 -0.92553991  
### 4.4 Construct cell type-specific gene network  
`constr_specific_net` constructs cell type-specific gene network from a given referenced gene interaction network by using the cell type scores obtained from `calc_score()`.

       net_specific = constr_specific_net(net = refnet,genescore = gene_score)
       
        > str(net_specific)  
            List of 8  
            $ Glutamatergic neuron          :List of 2  
             ..$ specific_net:'data.frame':	108942 obs. of  2 variables:  
             .. ..$ V1: chr [1:108942] "1000" "1000" "1000" "1000" ...  
             .. ..$ V2: chr [1:108942] "1001" "1004" "10076" "10092" ...  
             ..$ gene_score  :'data.frame':	4395 obs. of  2 variables:  
             .. ..$ gene : Factor w/ 4395 levels "1000","10000",..: 1637 2157 3203 984 703 2969 1183 1129 1128 103 ...  
             .. ..$ score: num [1:4395] 1.619 0.224 1.044 0.213 0.125 ...  
### 4.5 Identify cell type-specific disease gene module 
`constr_disease_module` identifies cell type-specific disease gene module from the cell type-specific gene network obtained from `constr_specific_net` and assesses its significance by permuation tests. Parameter 'returnALLmodules' is set to determine if all candidate cell type-specific disease gene modules should be returned or just return the significant ones.

      Cts_DGM = constr_disease_module(Cts_list = net_specific, disease_gene = disease_gene,
                                     padjCutoff = 0.1, returnAllmodules = FALSE, itera = 1000)
       
      > str(Cts_DGM)
      List of 8
       $ Glutamatergic neuron          :List of 1
        ..$ 1:List of 4
        .. ..$ disease_module:'data.frame':	1717 obs. of  2 variables:
        .. .. ..$ V1: chr [1:1717] "1006" "1006" "1006" "1006" ...
        .. .. ..$ V2: chr [1:1717] "11122" "1496" "1740" "1762" ...
        .. ..$ gene_score    :'data.frame':	238 obs. of  2 variables:
        .. .. ..$ gene : Factor w/ 4395 levels "1000","10000",..: 103 2217 3683 696 1110 413 2436 920 1549 1551 ...
        .. .. ..$ score: num [1:238] 0.6235 0.9655 0.0388 0.1092 0.317 ...
        .. ..$ Pvalue        : num 0
        .. ..$ p.adjust      : num 0
       $ GABAergic interneuron         :List of 1
        ..$ 1:List of 4
        .. ..$ disease_module:'data.frame':	2201 obs. of  2 variables:
        .. .. ..$ V1: chr [1:2201] "10055" "10055" "1006" "1006" ...
        .. .. ..$ V2: chr [1:2201] "158" "25942" "11122" "1995" ...
        .. ..$ gene_score    :'data.frame':	305 obs. of  2 variables:
        .. .. ..$ gene : Factor w/ 6395 levels "10001","10003",..: 142 3079 340 747 970 1545 594 1908 3380 65 ...
        .. .. ..$ score: num [1:305] 0.604 1.069 0.118 0.18 0.28 ...
        .. ..$ Pvalue        : num 0
        .. ..$ p.adjust      : num 0
### 4.6 (Optional) Export cell type-specific disease gene module to cytoscape
Export cell type-specific disease gene module to cytoscape for visualization, take `Glutamatergic neuron` for instance:

     exportNetToCytoscape<-Cts_DGM[["Glutamatergic neuron"]][[1]]$disease_module  
     table.score <- Cts_DGM[["Glutamatergic neuron"]][[1]]$gene_score  
     write.csv(exportNetToCytoscape,"disease_Module.csv",row.names = FALSE,quote = FALSE,row.names = FALSE,sep = '\t')  
     write.csv(table.score,"genescore.csv",row.names = FALSE,quote = FALSE,row.names = FALSE,sep = '\t')  
     
    > head(exportNetToCytoscape)  
              V1    V2  
    26629311 1006 11122  
    26629581 1006  1496  
    26629719 1006  1740  
    26629727 1006  1762  
    26629779 1006  1995  
    26630013 1006  2334  
    > head(table.score)  
         gene      score  
    10  10347 0.62349634  
    38  51412 0.96554854  
    41     81 0.03876118  
    87    158 0.10924185  
    93   2334 0.31703320  
    95 116986 1.57486599  


