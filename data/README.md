# Dataset Information: TCGA-LUAD PanCanAtlas

This project analyzes the **Lung Adenocarcinoma (LUAD)** cohort from **The Cancer Genome Atlas (TCGA) Pan-Cancer Atlas 2018** study.

## Preprocessed Dataset Included

For immediate exploration and modeling, the repository includes:
* `data/integrated_luad_dataset.csv` (78.6 KB, 417 patients, 21 multi-omics features)

This fused table contains clinical phenotypes, binary somatic driver mutation flags, total mutation burdens, and log2-transformed mRNA-seq RSEM expression levels.

---

## Raw Dataset Acquisition (For Reproducing ETL from Scratch)

To re-run the raw ETL ingestion pipeline from scratch (`src/etl.py`):
1. Visit the [cBioPortal Data Hub](https://www.cbioportal.org/datasets) or download directly from [cBioPortal LUAD PanCanAtlas](https://cbioportal-datahub.s3.amazonaws.com/luad_tcga_pan_can_atlas_2018.tar.gz).
2. Extract the archive into `data/luad_tcga_pan_can_atlas_2018/`.
3. Ensure the following files are present:
   * `data_clinical_patient.txt`
   * `data_mutations.txt`
   * `data_mrna_seq_v2_rsem.txt`
4. Run:
   ```bash
   python src/etl.py --raw_dir data/luad_tcga_pan_can_atlas_2018 --output data/integrated_luad_dataset.csv
   ```

*Note: Raw cBioPortal files (>350 MB) are excluded from this Git repository via `.gitignore` to prevent repository bloat.*
