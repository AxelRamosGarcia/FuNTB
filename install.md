# Installation of FunTB V.0.0.1

Last updated: Jan 31, 2024

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
git clone https://github.com/ind-genomics/FUN-TB.git
export FunTB_DIR=$(realpath FunTB/)
# You can put this export command in the your .bashrc file
# so that you don't need to type every time you run the FunTB
```

## Test run

Run FunTB scripts with test dataset

#### 1. Variation dictionary creation script

Once located within the FunTB directory the first step is to generate the variation dictionary file from MTBSeq V.0.1 output, in order to do this, run the following command:
```bash
python FunTB_dictionary.py MTBSeq_HIV_Data.tab Clinical_Data_HIV.csv
```
After the execution of this script a TXT file will be generated in the *`Variations_dictionaries`* which will contain the information of every sample and those genes that present any alteration.

#### 2. Phenotype-based samples lists generation script

Then to generate the samples' lists, run the following command:
```bash
python Sample_Grouping_Creation.py Clinical_Data_HIV.csv
```
The execution of this script will generate a series of TXTs files which will contain the samples of each of the groups with common shared clinical features.

#### 3. Phenotype-centric and gene-surrounded networks structuration script

Finally, to generate the network files, run the following command:
```bash
python FunTB.py VIH_network MTBSeq_HIV_Data.txt 30 0.2 0.1 0.7 VIH_positive_Samples.txt VIH_negative_Samples.txt
```
This will generate three XML-network format files which you can find within the Networks_files directory.

The FunTB full analysis took about 5 minutes by Intel(R) Xeon(R) CPU E5-2676 v3 @ 2.40GHz with 8 CPU cores.
