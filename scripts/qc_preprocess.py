import scanpy as sc
import os

DATA_DIR = "data"

def main():
    mtx_path = os.path.join(DATA_DIR, "matrix.mtx")
    if not os.path.exists(mtx_path):
        raise FileNotFoundError(f"{mtx_path} not found. Synthetic data should be in the data/ directory.")

    # Read synthetic 10X-style matrix
    adata = sc.read_10x_mtx(DATA_DIR, var_names="gene_symbols", cache=True)

    # Basic QC
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    adata = adata[
        (adata.obs["n_genes_by_counts"] > 1)
        & (adata.obs["n_genes_by_counts"] < 1000)
    ].copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0, max_mean=10, min_disp=0.0)
    if adata.var["highly_variable"].sum() > 0:
        adata = adata[:, adata.var["highly_variable"]].copy()

    os.makedirs("results", exist_ok=True)
    adata.write("results/preprocessed_anndata.h5ad")
    print("Saved preprocessed data to results/preprocessed_anndata.h5ad")

if __name__ == "__main__":
    main()
