"""
TCGA-LUAD Multi-Omics Data Ingestion and ETL Pipeline.

Extracts, cleans, standardizes, and fuses clinical, somatic mutation,
and mRNA-seq gene expression matrices from The Cancer Genome Atlas (TCGA)
Lung Adenocarcinoma (LUAD) PanCanAtlas cohort.
"""

import os
import argparse
import numpy as np
import pandas as pd

DRIVER_GENES = ["EGFR", "KRAS", "TP53", "STK11", "KEAP1"]


def clean_clinical_data(clinical_path: str) -> pd.DataFrame:
    """Load and standardize clinical patient records."""
    df_raw = pd.read_csv(clinical_path, sep="\t", comment="#")
    clinical = df_raw[[
        "PATIENT_ID",
        "AGE",
        "SEX",
        "AJCC_PATHOLOGIC_TUMOR_STAGE",
        "PATH_T_STAGE",
        "PATH_N_STAGE",
        "PATH_M_STAGE",
        "OS_MONTHS",
        "OS_STATUS",
        "DFS_STATUS"
    ]].copy()

    clinical["SAMPLE_ID"] = clinical["PATIENT_ID"] + "-01"
    clinical = clinical.rename(columns={
        "AGE": "Diagnosis Age",
        "AJCC_PATHOLOGIC_TUMOR_STAGE": "AJCC Pathologic Stage",
        "PATH_T_STAGE": "AJCC Pathologic T-Stage",
        "PATH_N_STAGE": "AJCC Pathologic N-Stage",
        "PATH_M_STAGE": "AJCC Pathologic M-Stage",
        "OS_MONTHS": "Overall Survival (Months)",
        "OS_STATUS": "Overall Survival Status",
        "DFS_STATUS": "Disease Free Status",
        "SEX": "SEX"
    })
    return clinical.drop(columns=["PATIENT_ID"])


def process_mutations(mutation_path: str, driver_genes=DRIVER_GENES) -> pd.DataFrame:
    """Extract binary somatic driver mutation matrix and total mutation counts."""
    mutation_raw = pd.read_csv(mutation_path, sep="\t", comment="#")
    mutation = mutation_raw[["Hugo_Symbol", "Tumor_Sample_Barcode"]].copy()
    mutation = mutation.rename(columns={"Tumor_Sample_Barcode": "SAMPLE_ID"})

    mutation_driver = mutation[mutation["Hugo_Symbol"].isin(driver_genes)]
    mutation_matrix = (
        mutation_driver.assign(mutated=1)
        .pivot_table(
            index="SAMPLE_ID",
            columns="Hugo_Symbol",
            values="mutated",
            aggfunc="max",
            fill_value=0
        )
        .reset_index()
    )
    mutation_matrix.columns = ["SAMPLE_ID"] + [
        f"{c}_mutated" for c in mutation_matrix.columns[1:]
    ]

    mutation_burden = mutation.groupby("SAMPLE_ID").size().reset_index(
        name="total_mutation_count"
    )

    return mutation_matrix.merge(mutation_burden, on="SAMPLE_ID", how="left")


def process_expression(expression_path: str, driver_genes=DRIVER_GENES) -> pd.DataFrame:
    """Extract, filter, and log2-normalize mRNA-seq RSEM expression values."""
    expr_raw = pd.read_csv(expression_path, sep="\t", comment="#")
    expr = expr_raw.dropna(subset=["Hugo_Symbol"]).copy()
    expr = expr[expr["Hugo_Symbol"].isin(driver_genes)]
    expr = expr.drop_duplicates(subset=["Hugo_Symbol"])
    expr = expr.set_index("Hugo_Symbol")

    if "Entrez_Gene_Id" in expr.columns:
        expr = expr.drop(columns=["Entrez_Gene_Id"])

    expr = expr.transpose()
    expr.index.name = "SAMPLE_ID"
    expr = expr.reset_index()
    expr.columns = ["SAMPLE_ID"] + [f"{g}_expr" for g in expr.columns[1:]]

    for c in expr.columns[1:]:
        expr[c] = pd.to_numeric(expr[c], errors="coerce")
        expr[c] = np.log2(expr[c] + 1)

    return expr


def run_pipeline(raw_dir: str, output_path: str) -> pd.DataFrame:
    """Execute end-to-end multi-omics ETL pipeline."""
    clinical_file = os.path.join(raw_dir, "data_clinical_patient.txt")
    mutation_file = os.path.join(raw_dir, "data_mutations.txt")
    expression_file = os.path.join(raw_dir, "data_mrna_seq_v2_rsem.txt")

    print("[*] Processing clinical data...")
    clinical = clean_clinical_data(clinical_file)

    print("[*] Processing somatic mutations...")
    mutation = process_mutations(mutation_file)

    print("[*] Processing mRNA-seq expression...")
    expr = process_expression(expression_file)

    print("[*] Fusing multi-omics modalities...")
    integrated = clinical.merge(mutation, on="SAMPLE_ID", how="inner") \
                         .merge(expr, on="SAMPLE_ID", how="inner")

    integrated = integrated.drop_duplicates(subset="SAMPLE_ID")

    print("[*] Imputing missing values...")
    for col in integrated.columns:
        if integrated[col].dtype in ["float64", "int64"]:
            integrated[col] = integrated[col].fillna(integrated[col].median())
        else:
            integrated[col] = integrated[col].fillna("Unknown")

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    integrated.to_csv(output_path, index=False)
    print(f"[✓] Pipeline complete! Integrated dataset saved to: {output_path}")
    print(f"    Dimensions: {integrated.shape[0]} patients x {integrated.shape[1]} features")
    return integrated


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TCGA-LUAD Multi-Omics ETL Pipeline")
    parser.add_argument("--raw_dir", type=str, default="data/luad_tcga_pan_can_atlas_2018",
                        help="Directory containing unzipped cBioPortal TCGA-LUAD files")
    parser.add_argument("--output", type=str, default="data/integrated_luad_dataset.csv",
                        help="Path to save integrated dataset CSV")
    args = parser.parse_args()

    run_pipeline(args.raw_dir, args.output)
