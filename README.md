# FuNTB

*`FuNTB`* is a free, open-source Python tool designed for analyzing output from MTBSeq v.1.0. This standalone command-line tool enables the comparison of Single Nucleotide Polymorphisms (SNPs) across phenotypically diverse sets of Mycobacterium tuberculosis samples. By structuring the data into phenotype-centered networks, FUN-TB allows users to visualize which genes are altered and how these alterations correlate with specific phenotypes. The tool provides insights into gene alterations by representing them as nodes, with node size and edge width reflecting the strength of their relationship to particular phenotypes. 

*`FuNTB`* consists of three core scripts:

1. **Variation Dictionary Creation**: Generates a dictionary of variations found in the `MTBSeq` output.
2. **Phenotype-Based Sample List Generation**: Generates a dictionary of variations found in the `MTBSeq` output.
3. **Phenotype-Centric and Gene-surrounded Network Structuring**: Constructs networks based on the relationship between phenotypic traits and gene alterations, making it easy to identify key genetic markers associated with specific phenotypes.

## Variation dictionary creation script
This script aims to create a Python dictionary object from an MTBSeq V.1.0 file. The dictionary summarizes the presence of SNPs, along with the position and frequency parameters of each altered gene per sample.



## Phenotype-based samples lists generation script
This second script aims to generate a series of lists containing sample IDs. These lists are constructed such that samples within the same list share common phenotypical features and clinical values. The input CSV file must follow the format below to ensure that the clinical data is processed correctly.

### Required Columns:
- **Sample_ID**: A unique identifier for each sample.

### Optional Metadata Columns (Can vary depending on the dataset but should be consistent within a file):
- **Age:** Age of the individual (numeric).
- **Gender:** Gender of the individual (e.g., `Male`, `Female`, `Other`).
- **Year_of_sample:** Year when the sample was collected (e.g., `2020`, `2021`, `2022`, `2023`).
- **Sequencing_technology:** The sequencing technology used (e.g., `Illumina`, `Nanopore`, `PacBio`).
- **Geographic_location:** Location where the sample was collected (e.g., `Mexico`, `USA`, `Spain`, `Germany`).
- **Disease_status:** Clinical status of the individual (e.g., `HIV Positive`, `Tuberculosis`, `Healthy`).
- **Treatment:** Treatment the individual is receiving (if applicable, e.g., `Antibiotic`, `Antiretrovirals`, `None`).
- **Other_metadata:** Any additional information relevant to the study.

#### Example Clinical Data File (`sample_clinical_data.csv`):
| Sample_ID | AGE | GENDER | YEAR_OF_SAMPLE | SEQUENCING_TECHNOLOGY | GEOGRAPHIC_LOCATION | DISEASE_STATUS  | TREATMENT       |
|-----------|-----|--------|----------------|-----------------------|---------------------|-----------------|----------------|
| S1        | 45  | Male   | 2021           | Illumina              | USA                 | HIV positive    | Antiretrovirals |
| S2        | 60  | Female | 2019           | Nanopore              | Brazil              | Tuberculosis    | Antibiotics     |
| S3        | 30  | Male   | 2020           | PacBio                | Germany             | Healthy         | None           |
| S4        | 50  | Female | 2022           | Illumina              | Canada              | HIV positive    | Antiretrovirals |

#### File Requirements:
1. The file must be in CSV format.
2. The first column (`Sample_ID`) is mandatory and should contain unique identifiers.
3. Other columns can vary based on the dataset, but all values within a column should be consistent (e.g., no mixing of strings and numeric values in the same column).
4. Ensure there are no extra spaces or typos in column names or values.

#### Flexible processing
The script dynamically adapts to different metadata columns. Any metadata category can be used for grouping and filtering if the column name is correctly provided when running the script.

## Phenotype-centric and gene-surrounded networks structuration script.

## Table of Contents

- [Features](#features)
- [Set Up](#Setup)
- [Usage](#usage)
- [Documentation](#documentation)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)
