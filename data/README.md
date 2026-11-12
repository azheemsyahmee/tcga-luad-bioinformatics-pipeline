# Dataset Information: TCGA-LUAD PanCanAtlas

This project analyzes the **Lung Adenocarcinoma (LUAD)** cohort from **The Cancer Genome Atlas (TCGA) Pan-Cancer Atlas 2018** study, distributed via the cBioPortal Data Hub.

---

## 1. Data Access Tier & Open-Access Scope

Under the **National Cancer Institute (NCI) Genomic Data Commons (GDC) Data Access Policy**, TCGA datasets are strictly categorized into two tiers:

* **Open Access Tier (Used in this repository):**
  - Includes de-identified clinical metadata, curated somatic mutation annotations (MAF files without germline identifiers), and normalized gene-level mRNA-seq expression values (RSEM quantification).
  - Open-access data are publicly available for research and education without requiring dbGaP approval.
  - **All data ingested, processed, and bundled in this project belong strictly to this Open Access tier.**
* **Controlled Access Tier (Not used or included):**
  - Includes individual-level germline variant calls, raw sequencing reads (`.bam` / `.fastq`), high-density SNP arrays, and sensitive clinical identifiers that carry re-identification risk.
  - Controlled-access data require an approved application through NIH dbGaP and NCI Data Access Committee (DAC) authorization.
  - **This repository does not use, store, or redistribute any Controlled Access TCGA data.**

---

## 2. Preprocessed Dataset Included

For immediate exploration and modeling, the repository bundles:
* `data/integrated_luad_dataset.csv` (78.6 KB, 417 patients, 21 multi-omics features)

This fused table contains de-identified clinical phenotypes, binary somatic driver mutation flags, total mutation burdens, and log2-transformed mRNA-seq RSEM expression levels for major lung adenocarcinoma driver genes.

---

## 3. Raw Source Data & ETL Acquisition

The raw input tables originate from the **cBioPortal Pan-Cancer Atlas 2018** archive (`luad_tcga_pan_can_atlas_2018`):
1. **Clinical Data:** `data_clinical_patient.txt` — de-identified patient demographic and survival records.
2. **Somatic Mutations:** `data_mutations.txt` — public somatic variant calls annotated by cBioPortal.
3. **Expression Profiles:** `data_mrna_seq_v2_rsem.txt` — normalized gene-level mRNA-seq RSEM abundance estimates.

To re-run the raw ETL ingestion pipeline from scratch (`src/etl.py`):
1. Download the archive from the [cBioPortal Data Hub](https://www.cbioportal.org/datasets) or directly from [cBioPortal LUAD PanCanAtlas](https://cbioportal-datahub.s3.amazonaws.com/luad_tcga_pan_can_atlas_2018.tar.gz).
2. Extract the archive into `data/luad_tcga_pan_can_atlas_2018/`.
3. Verify the three files above are present.
4. Execute the pipeline:
   ```bash
   python src/etl.py --raw_dir data/luad_tcga_pan_can_atlas_2018 --output data/integrated_luad_dataset.csv
   ```

*Note: Raw cBioPortal archives (>350 MB) are excluded from this Git repository via `.gitignore`.*
