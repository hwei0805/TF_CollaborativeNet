################# Installation  ###########

1. Save the package on the server.

2. extract the package by typing:

tar zxvf SCCM_v1.1.tar.gz 

3. install the modules by tying the following commands:

cd SCCM_v1.1

sudo install.sh

Note: You need to have the admin privileges to install the modules. This will install  Config::Simple,Statistics::RankCorrelation and Parallel::ForkManager total 3 perl cpan modules.

If it does not work, cd to individual module directories and follow the directions to install the three modules individually. 

4. Edit SCCM_pipe.cfg: replace the lib=pathToscripts to the full directory path to scripts directory.

5. add the scripts directory to your path.

############ how to run ########

1. make a copy of SCCM_pipe.cfg, you could copy it to the directory where your two input files are and modify it. Usually you only need to change the number of CPUs and the two input files names.

2. type:

nohup perl SCCM_pipe.pl SCCM_pipe.cfg & 

########### about config file ##############

Several pipeline parameters can be set through config file

1. cpu:  how many cpus you want to use in the analysis. The value should be equal or less than the total cpus your machine has.

2. topPick: number of top genes you want to examine for the building of the connectivity matrix. Refer to the paper for details. Default 100. Reducing the number will result in smaller clusters.

3. tripleLinkN(N=1,2,3), again refer to the paper for a description of these three parameters.  Higher values will result in smaller clusters. The values should be set such as: tripleLink1 > tripleLink2 > tripleLink3. 
 

######## input files format #########

Two input files 

1. expression file, title line should start with # eg  (C1--chip data set 1, RNA-seq data set is fine too)

#gene  C1   C2   C3   C4   

AAM    2.98     3.45   5.33   3.22

2 Transcription factor (TF)  gene list, one gene per line. The TF gene identifier should be the same as used in the gene expression file. All TF genes should be included in the gene expression file. 

3. Note: The pipeline will detect the matrix file. If it is exists, the  pipeline will skip the co-expression matrix building step. Remove the matrix file before re-runnning the pipeline if the existing matrix is not what you want. 

4. SCCM_builder.pl and SCCM_decomposition.pl could be run separately, type program without any parameter for help info. 


########## output files ###############

1.The final cluster file name is started with Cluster_. Clusters are sorted by maximum connectivity. 

2.  Intermediate  files:

2.1 .coexp.matrix.txt file. The row names are TFgene_1,TFgene_2,TFgene_3 etc....  TFgene_n, the col names are also TFgene_1, TFgene_2, TFgene_3, etc.... The rows and columns are in the same order. The value is the number of top 100 shared genes for TFgene_x and TFgene_y. The maximum value is 100 and minimum is 0.

2.2 top_100 directory includes the top 100 co-expressed genes for each TFgene.  The file name is the TFgene identifier.

2.3 The detailed TFgenes connectivity in each cluster is included in the cluster_within_connectivity directory. The format is TFgene1   TFgene2  connectivity(number of top 100 gene shared). This file can be used as input to Cytoscape or other tools to generate the connectivity graph.

######### help info #############

For help contact Jeff Nie at jnie@morgridgeinstitute.org or 

                 Hairong Wei at hairong@mtu.edu


