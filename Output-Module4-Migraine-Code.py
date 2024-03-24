# ORGINAL MAPPING!

# CHR -> ['chromosome']
# SNP -> ['variant_id']
# BP -> ['base_pair_location', 'genpos']
# A1 -> ['effect_allele']
# A2 -> ['other_allele']
# N -> []
# P -> ['p_value']
# INFO -> ['info']
# MAF -> ['effect_allele_frequency']
# BETA -> ['beta', 'p_bolt_lmm_inf']
# SE -> ['standard_error']
# OR -> []
# Z -> []
# D -> []
# Transformation for code !

# chromosome -> CHR
# variant_id -> SNP
# base_pair_location -> BP
# genpos -> BP
# effect_allele -> A1
# other_allele -> A2
# NA_N -> N
# p_value -> P
# info -> INFO
# effect_allele_frequency -> MAF
# beta -> BETA
# p_bolt_lmm_inf -> BETA
# standard_error -> SE
# NA_OR -> OR
# NA_Z -> Z
# NA_D -> D
 Here's how you can perform the required conversion using pandas in Python:
```python
import pandas as pd

# Load the original GWAS data into a pandas DataFrame
df = pd.read_csv('gwas.csv')

# Perform the required column renaming and reordering
df_renamed = df.rename(columns={'chromosome': 'CHR', 'variant_id': 'SNP', 'base_pair_location': 'BP',
                               'genpos': 'BP', 'effect_allele': 'A1', 'other_allele': 'A2',
                               'p_value': 'P', 'info': 'INFO', 'effect_allele_frequency': 'MAF',
                               'beta': 'BETA', 'p_bolt_lmm_inf': 'BETA', 'standard_error': 'SE'})
df_reordered = df_renamed[['CHR', 'SNP', 'BP', 'A1', 'A2', 'P', 'INFO', 'MAF', 'BETA', 'SE']]

# Save the modified DataFrame to a new CSV file
df_reordered.to_csv('finalgwas.csv', index=False)
```
This code will load the `gwas.csv` file into a pandas DataFrame, rename the specified columns, reorder them, and then save the resulting DataFrame to a new file called `finalgwas.csv`. The `index=False` argument ensures that no row indices are saved to the new file.