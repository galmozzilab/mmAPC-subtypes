this folder contains regulon scores used to merge with the seurat object. 
this is the final output from the scenic pipeline ran from their docker container.

scenic run steps for BAT and WAT:

#1. generate adjacencies matrix- grnboost 
expr_matrix = pd.read_csv('regressed/normalized_counts.csv', index_col=0) #this is where we saved our matrix extracted from seurat

### Transpose the matrix
expr_matrix_transposed = expr_matrix.transpose()

### Save the transposed matrix to a TSV file
expr_matrix_transposed.to_csv('expr_mat.tsv', sep='\t')

### load in motifs
MOTIFS_MGI_FNAME = 'motifs-v9-nr.mgi-m0.001-o0.0.tbl' #this one we downloaded from aertslab.org
OUT_TFS_MGI_FNAME = 'allTFs_mm.txt'

df_motifs_mgi = pd.read_csv(MOTIFS_MGI_FNAME, sep='\t')
mm_tfs = df_motifs_mgi.gene_name.unique()
with open(OUT_TFS_MGI_FNAME, 'wt') as f:
    f.write('\n'.join(mm_tfs) + '\n')
len(mm_tfs)

docker run -it --rm `
    -v "F:\scenic\bat\data:/data" `
    aertslab/pyscenic:0.12.1 pyscenic grn `
        --num_workers 4 `
        -o /data/expr_mat.adjacencies.tsv ` #this is the output file
        /data/expr_mat.tsv `
        /data/allTFs_mm.txt

# 2. run ctx step, will get regulons_v10.csv
docker run -it --rm `
    -v "F:\scenic\wat\data:/data" `
    aertslab/pyscenic:0.12.0 pyscenic ctx `
        /data/expr_mat.adjacencies.tsv `
        /data/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather ` #this one we downloaded from arertslab.org
        /data/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather ` #this one we downloaded from arertslab.org
        --annotations_fname /data/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl `
        --expression_mtx_fname /data/expr_mat.tsv `
        --mode "custom_multiprocessing" `
        --output /data/regulons_v10.csv `
        --num_workers 4

#3. run aucell module after generated prerequisite files (regulons_v10)
### will generate auc_mtx_v10.csv
docker run -it --rm `
    -v "F:\scenic\wat\data:/data" `
    aertslab/pyscenic:0.12.0 pyscenic aucell `
        /data/expr_mat.tsv `
        /data/regulons_v10.csv `
        -o /data/auc_mtx_v10.csv `
        --num_workers 4