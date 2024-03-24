# GWASPokerforPRS

GWAS Summary Statistic Tool: A Meta-Analysis and Parsing Tool for Polygenic Risk Score Calculation

GWAS (genome-wide association study) summary statistic files are used to calculate polygenic risk scores (PRS). Multiple research groups provide these files for a specific phenotype or disease. Scanning the GWAS available in the GWAS Catalog revealed that for a particular disease, there can be multiple GWAS files with varying populations, numbers of associations identified, sample sizes, validation samples, types of analyses used to generate the GWAS, genome builds, and types of information listed in the GWAS. Finding, downloading, and verifying the GWAS file for a specific phenotype can be challenging, as it involves downloading, extracting, parsing, and manually scanning the files. We propose a tool generated after analyzing 60,400 GWAS summary files from the GWAS Catalog. It allows scanning the GWAS file without fully downloading it. The process involves searching and downloading metadata for all GWAS related to a specific phenotype, partially downloading the GWAS file, parsing and cleaning the file to create a readable table, searching for specific column headers, and listing the necessary columns for polygenic risk score calculation. It also includes extracting the DOI and citation of the article using PMID, along with a Python code generator module for mapping original GWAS columns to those required by PRS tools.

### Module 1
```bash
python Module1-SearchPhenotypeandPopulation.py --phenotype asthma --population European
python Module1-SearchPhenotypeandPopulation.py --phenotype migraine

### Module 2
```bash
Once you manually process the output file, do not change the file name and keep it the same for smooth working of the code.
python Module2-Search_Poke_Normalize_Scan.py --processedfile migraine.csv

### Module 3
```bash
python Module3-DownloadGWAS.py --processedfile Input-Module3-Migraine.csv --indexer 1

The Input-Module3-Migraine.csv contains the following information. The name is the directory in which the file should be downloaded and the further processed files will be stored in the same directory.
Name,Download Link
migraine,http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90038001-GCST90039000/GCST90038646/GCST90038646_buildGRCh37.tsv

### Module 4
```bash
python Module4-ExtractGWAS.py --processedfile Input-Module3-Migraine.csv --indexer 1

### Module 5
```bash
If you want to check the columns in your own GWAS file. The code accepts .csv format.
python Module5-ListPRSColumns.py --gwasfile gwas.csv.modified


