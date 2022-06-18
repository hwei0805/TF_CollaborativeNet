# TF_CollaborativeNet (Previous called TF-Cluster): 

Note that this pipeline contain the most efficient software pipeline for identifying transcription factors (TFs) that work together do a job. 

The pipeline (in Perl) uses the high-throughput gene expression data (from the same tissue, RNA-seq or microarray) to build a coilllaborative subnetworks, each of them controls a biological process or a complext trait.

The pipeline does it in two steps

Step 1: Building a global collabroative network of all transcription factors (TFs) through genomew-wide coexpression analysis. Based on coexpression analysis results, it constructs TF collabroative network;

Step 2: Decompose the global collabroative network into many subnetworks. The subnetworks output earlier have tighter collaboration that later one.  Usually top 30 subnetworks are more biologically meaningful and can be interpreted eassily.

When you get a subnetwork, study the genes in each cluster, you will find they are TFs that control the same biological pathway/process or a complex trait. 

Read the publications below for more details.


######  How to cite ? ########

1) Nie, J., R. Stewart, F. Ruan, J. Thomson, H. Zhang, X. Cui and H. Wei*. 2011. TF-Cluster: a pipeline for identifying functionally coordinated transcription factors via network decomposition of the shared coexpression connectivity matrix (SCCM). BMC Systems Biology, 5:53. https://doi.org/10.1186/1752-0509-5-53.

2)	Ji, X. S. Chen, J. Li, W. Deng, Z. Wei and H. Wei*.  2017.  SSGA and MSGA: two seed-growing algorithms for constructing collaborative subnetworks. Scientific Reports. 7:1446, DOI:10.1038/s41598-017-01556-z
