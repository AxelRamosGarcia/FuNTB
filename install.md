# Installation of FunTB V.0.0.1

Last updated: Feb 23, 2025

FunTB is an open source standalone tool available for academic purposes. For commercial use of FunTB please contact: industrial.genomics@gmail.com

Please feel free to post on **Issues** or contact axel.ramos3737@gmail.com

## 1. FunTB requirements

 * Python: 3.11.5
 * Pandas: 2.0.2
 * Numpy: 1.24.3
 * Matplotlib: 3.7.1
 * Networkx: 3.1
 * Seaborn: 0.13.0
 * Scipy: 1.11.1

### 1.1 Install dependencies

##### Python v.3.11.5

[Windows](https://www.python.org/ftp/python/3.11.5/python-3.11.5-amd64.exe)
[Linux](https://www.python.org/ftp/python/3.11.5/Python-3.11.5.tgz)
[macOS](https://www.python.org/ftp/python/3.11.5/python-3.11.5-macos11.pkg)

##### Pandas 2.0.2
```bash
# Pip: Run pip install pandas == <version>
pip install pandas == 2.0.2
```
```bash
# Anaconda: Run conda install pandas = <version>
conda install pandas = 2.0.2
```
##### Numpy 1.24.3
```bash 
# Pip: Run pip install numpy == <version>
pip install numpy == 1.24.3
```
```bash
# Anaconda: Run conda install numpy = <version>
conda install numpy = 1.24.3
```

> [!Note]
> If you have already installed a specific version of NumPy and want to install another version, use "--force-reinstall" e.g., "pip install -- force-reinstall numpy == version".
> 
> if you have a virtual env and want to override the installation without uninstalling, use "--ignore-installed", e.g., "pip install --ignore-installed numpy == version".

##### Matplotlib 3.7.1
```bash 
# Pip: Run pip install matplotlib == <version>
pip install matplotlib == 1.24.3
```
```bash
# Anaconda: Run conda install matplotlib = <version>
conda install matplotlib = 1.24.3
```
##### Networkx 3.1
```bash 
# Pip: Run pip install networkx == <version>
pip install networkx == 3.1
```
```bash
# Anaconda: Run conda install networkx = <version>
conda install networkx = 3.1
```
##### Seaborn 0.13.0
```bash 
# Pip: Run pip install seaborn == <version>
pip3 install seaborn
```
```bash
# Anaconda: Run conda install conda-forge::seaborn
conda install conda-forge::seaborn
```
##### Scipy 1.11.1
```bash 
# Pip: Run python -m  pip install scipy
python -m pip install scipy
```
```bash
# Anaconda: Run conda install scipy
conda install scipy
```
## 2. Download FunTB

Use GitHub clone to download

```bash
cd $HOME  # or wherever you want
git clone https://github.com/AxelRamosGarcia/FuNTB.git
export FunTB_DIR=$(realpath FunTB/)
# You can put this export command in the your .bashrc file
# so that you don't need to type every time you run the FunTB
```

## Test run

This dataset contains 462 clinical samples from individuals with tuberculosis (TB), annotated with key phenotypic and clinical characteristics relevant to drug resistance profiling. The dataset includes:

### Clinical covariates:

- HIV: HIV co-infection status
- Diabetes: Diabetes comorbidity status

### Drug resistance phenotypes:
- Binary indicators (0 = susceptible, 1 = resistant) for 14 anti-TB drugs:
- DLM (Delamanid), CFZ (Clofazimine), MXF (Moxifloxacin), BDQ (Bedaquiline), INH (Isoniazid), LZD (Linezolid), RIF (Rifampicin), RFB (Rifabutin), LEV (Levofloxacin), EMB (Ethambutol), ETH (Ethionamide), PAS (Para-aminosalicylic acid), KAN (Kanamycin), AMI (Amikacin).

Purpose
This dataset is designed to support analyses of:

- Associations between clinical factors (e.g., HIV, diabetes) and drug resistance.
- Co-resistance patterns across first- and second-line TB drugs.
- Predictive modeling of resistance to guide treatment regimens.

Format
- Structured as a CSV file with rows representing samples and columns as variables.
- Binary phenotypes (S/R).

### 1. Variation dictionary creation script

Overview:
This test dataset for Script 1 comprises a sample MTBseq output file (or an equivalent annotated VCF file) that contains pre-filtered, non-synonymous SNPs from a subset of clinical TB samples. The file is structured with columns indicating the sample ID, gene name, mutation position, and mutation frequency.

Input Format:

- A tab-delimited file containing variant data as generated by MTBseq.
- A CSV file containing clinical data.

Expected Output:

- A Python dictionary (text file) where keys are sample IDs and values are nested dictionaries mapping genes to a further dictionary of mutation positions and their corresponding frequencies.
- For example:
```python
{
    "sample_001": {"rpoB": {450: 0.8, 452: 0.2}, "katG": {315: 1.0}},
    "sample_002": {"rpoB": {450: 0.9}}
}
```

Once located within the FunTB directory the first step is to generate the variation dictionary file from MTBSeq V.0.1 output, in order to do this, run the following command:
```bash
python FunTB_dictionary.py MTBseq_file.tab dictionary_clinical_test_data.csv
```
or

```bash
python FunTB_dictionary_vcf.py
```

After the execution of this script a TXT file will be generated in the *`Variations_dictionaries`* which will contain the information of every sample and those genes that present any alteration.

### 2. Clinical Sample Grouping Test Dataset

Overview:
This test dataset for Script 2 is a clinical metadata CSV file consisting of key clinical and phenotypic variables for each sample, such as HIV co-infection status, diabetes status, and binary indicators for resistance to 14 anti-TB drugs. This data supports the grouping of samples based on user-defined filtering conditions.

Input Format:

- A CSV file where each row represents a clinical sample.
- Mandatory columns include Sample ID (unique identifier) alongside additional columns for variables such as HIV, Diabetes, and drug resistance phenotypes (e.g., INH, RIF, etc.).

Expected Output:

- One or more text files, each containing a list of Sample IDs that match a set of user-defined filtering conditions.
- For instance, if the user defines a group for rifampicin-resistant samples by filtering on Resistance Status = Resistant for Rifampicin and Species = Mycobacterium tuberculosis, a text file such as Rifampicin_Resistant.txt will be generated.

Usage:

Then to generate the samples' lists, run the following command:
```bash
python Sample_Grouping_Creation.py test_clinical_data.csv DIABETES
```
The execution of this script will generate a series of TXTs files which will contain the samples of each of the groups with common shared clinical features.

### 3. Network Construction Test Dataset

Overview:
The test dataset for Script 3 includes the outputs from Scripts 1 and 2. This script integrates the variation dictionary and the grouped sample ID lists to generate an XML-based network file that visualizes the associations between genes and phenotypic groups.

Input Format:

- The variation dictionary (text file) produced by Script 1.
- Grouped sample lists (text files) generated by Script 2.

Expected Output:

- An XML file representing the gene-phenotype network, where nodes correspond to genes (with sizes reflecting CAIS values) and phenotype ego-nodes.
- The network file should be compatible with visualization tools such as Cytoscape or Gephi, enabling users to explore mutation patterns and their association with drug resistance phenotypes.

Usage:

Finally, to generate the network files, run the following command:
```bash
python FuNTB.py test_network test_dictionary.txt 3 0.9385527090157504 0.0007787658410143285 0.060668525143235286 test_positive_Samples.txt test_negative_Samples.txt
```
This will generate three XML-network format files which you can find within the Networks_files directory.

## 3. Environment Setup (Optional)

You can quickly set up a virtual environment using the provided *`funtb.yml`* file to ensure that all the dependencies are installed automatically.


### 1. Create the environment:
```bash
conda env create -f funtb.yml
```

### 2. Activate the environment
```bash
conda activate funtb
```

This will create a dedicated environment with all necessary dependencies pre-installed. You can now proceed with the usage steps.
