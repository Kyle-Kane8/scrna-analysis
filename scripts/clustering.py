import scanpy as sc
import os

def main():
    infile = "results/preprocessed_anndata.h5ad"
    if not os.path.exists(infile):
        raise FileNotFoundError("Run qc_preprocess.py first.")
    adata = sc.read(infile)

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=10)
    sc.pp.neighbors(adata, n_neighbors=5, n_pcs=10)
    sc.tl.umap(adata)
    sc.tl.tsne(adata)
    sc.tl.leiden(adata, resolution=0.5)
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")

    adata.write("results/clustered_anndata.h5ad")
    print("Saved clustered data to results/clustered_anndata.h5ad")

if __name__ == "__main__":
    main()
