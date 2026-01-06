import scanpy as sc
import os

def main():
    infile = "results/clustered_anndata.h5ad"
    if not os.path.exists(infile):
        raise FileNotFoundError("Run clustering.py first.")
    adata = sc.read(infile)

    os.makedirs("results/figures", exist_ok=True)

    sc.pl.umap(adata, color=["leiden"], save="_leiden.png", show=False)
    example_genes = list(adata.var_names[:2])
    for g in example_genes:
        sc.pl.umap(adata, color=g, save=f"_{g}.png", show=False)

    sc.pl.rank_genes_groups(adata, n_genes=5, sharey=False, save="_markers.png", show=False)
    print("Saved figures (Scanpy by default writes to `figures/`).")

if __name__ == "__main__":
    main()
