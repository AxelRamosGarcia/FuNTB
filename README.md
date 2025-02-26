# FuNTB

*`FuNTB`* is a free, open-source Python tool designed for analyzing output from MTBSeq v.1.0. This standalone command-line tool enables the comparison of Single Nucleotide Polymorphisms (SNPs) across phenotypically diverse sets of Mycobacterium tuberculosis samples. By structuring the data into phenotype-centered networks, FUN-TB allows users to visualize which genes are altered and how these alterations correlate with specific phenotypes. The tool provides insights into gene alterations by representing them as nodes, with node size and edge width reflecting the strength of their relationship to particular phenotypes. 

*`FuNTB`* consists of three core scripts:

1. *Variation Dictionary Creation*: Generates a dictionary of variations found in the `MTBSeq` output.
2. *Phenotype-Based Sample List Generation*: Generates a dictionary of variations found in the `MTBSeq` output.
3. *Phenotype-Centric and Gene-surrounded Network Structuring*: Constructs networks based on the relationship between phenotypic traits and gene alterations, making it easy to identify key genetic markers associated with specific phenotypes.


## Phenotype-based samples lists generation script
This second script aims to generate a series of lists containing sample IDs. These lists are constructed such that samples within the same list share common phenotypical features and clinical values. The input CSV file must follow the format below to ensure that the clinical data is processed correctly.

### Required Columns:
- *Sample_ID*: A unique identifier for each sample.

### Optional Metadata Columns (Can vary depending on the dataset but should be consistent within a file):
- *Age:* Age of the individual (numeric).
- *Gender:* Gender of the individual (e.g., `Male`, `Female`, `Other`).
- *Year_of_sample:* Year when the sample was collected (e.g., `2020`, `2021`, `2022`, `2023`).
- *Sequencing_technology:* The sequencing technology used (e.g., `Illumina`, `Nanopore`, `PacBio`).
- *Geographic_location:* Location where the sample was collected (e.g., `Mexico`, `USA`, `Spain`, `Germany`).
- *Disease_status:* Clinical status of the individual (e.g., `HIV Positive`, `Tuberculosis`, `Healthy`).
- *Treatment:* Treatment the individual is receiving (if applicable, e.g., `Antibiotic`, `Antiretrovirals`, `None`).
- *Other_metadata:* Any additional information relevant to the study.

