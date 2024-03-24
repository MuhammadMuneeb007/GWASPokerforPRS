import argparse
import pandas as pd
import numpy as np
import pandas as pd
import nltk
import re
import os
import pandas as pd
from fuzzywuzzy import fuzz
import requests
from bs4 import BeautifulSoup
import os
import shutil
import pandas as pd
import chardet

import pandas as pd
import chardet
import os
import subprocess
chromosome_list = ['chr:pos','chr:position','ch','var_name','name'
'chromosome_position_reference_allele_other_allele_b37','chromosome','chr',
'chr_build36','chr_pos_b36','chro','chrom','chromo']

snp_list = ['snp','markername','marker','rs','rsid','snpid','var_name',
'rs_number','dbsnprsid/marker','snp_ids','snp_id','snpname','name',
'rs_numbers','variant_id','rs_id','var_id','rs_marker_id',
'snprsid','illuminasnp','rs_id_all','rsnumber','rsq',
'marker_id','locus_id','variant_id_when_present','id','varid','marker_id',
'uniqid','oldid','dbsnprsid/marker','markerid','variantid','vcf_ids','pmarkername']

base_pair_list = ['name','var_name','name',
'base_pair_location','ase_pair_location','base_pair_locatio','base_pair_loca','bp','markername','marker_id','markerid','marker_id',
'pos','position','id','chr_pos_b36','pos_b37', 'phys_pos', 'positionhg19', 'position_build36', 'bpos', 'bpos', 'position_b37', 'posb37', 
'position_b37', 'position_grch37', 'pos_build36', 'chr:position', 
'base_pair_position', 'genpos', 'pos_b37', 'pos_b37', 'pos_b37', 
'position_hg19', 'position_hg19','chromosome_position_reference_allele_other_allele_b37',
'snp_pos','base_pair_location_grch37', 'e_pair_location', 'e_pair_location', 'pair', 'pair', 'base_pair_lo', 'base_pair_locations', 'base_pair', 'base_pair_location_grch37', 'base_pa', 'base_pa',
'base', 'base', 'base_pair_lo', 'base_', 'base_pair_locations', 'base_pair','locus_id','b'
]


effect_allele_list = [
'effect_allele','a1','allele1',
'allele_1','effect_allele','reference_allele','inc_allele','ea',
'a_1','allelea','ref','effect_allel','eff_allele','effect_all','minor_alele',
'alleles','effect_allele_plink','effect_allele_saige','minor_allele','t_allele',
'effect-allele','minorallelefrequency','effectallele','referenceallelesa/b','ref_allele',
'effectallele','allele.1','e','ea','allele_0','allele0','coded_allele','effect_allele_all'
'minorallele','effectallele','minorallele','a0'
]


alternative_allele_list = [
'a_2','a2', 'allele2', 'allele_2', 'other_allele', 'non_effect_allele', 'dec_allele', 'nea'
'alt_allele','alternate_allele','alternative_allele','alternate_ids','altfreq1','alt',
'alternative_allele','alt_allele','other-allele','otherallele','a2_other','other','other_all',
'other.all','allele1','alleleb','noneff_allele','allele_1','allele1','noneffect_allele','noncoded_allele'
'noneffect_allele','non_effect_allele_all','non_coded_allele','noneffectallele','majorallele',
'noneff_allele','oa','a1','none_effect_alllele','nea'
]

N_list = [
'n', 'ncase', 'cases_n', 'n_case', 'n_cases', 'n_controls', 'n_cas', 
'n_con', 'n_case', 'ncontrol', 'controls_n', 'n_control', 'weight',
'n_total','all_total', 'cases_total', 'chr_x_all_total', 'chr_x_all_total_females',
'chr_x_all_total_males', 'controls_total', 'n_total_sum', 'ntotal', 'total_n', 
'total_sample_size', 'total_samplesize','totalsamplesize',
'samplesize', 'n_samples', 'n_samples', 'mantra_n_samples', 'sample_size',
'n_effective_samplesize', 'n_effective_samplesize_all', 'n_effective_samplesize_females', 
'n_effective_samplesize_males', 'n_samples', 'samplesize', 'n_samples',
'samplesize', 'n_samples', 'sample_size', 'sample_size', 'sample_size',
'sample_n', 'effective_sample_size', 'sample-size', 'sample-size',
'sample-size-cases', 'n_samples', 'sample_size', 'samplesize', 'n_samples', 
'samplesize', 'n_samples', 'sample_num', 'n_samples', 'sample_size', 'nsample',
'n_samples', 'n_samples', 'n_samples', 'nsample', 'nsample', 'n_samples', 
'n_effective_samples_study_level', 'n_effective_samples_variant_level',
'n_samples_study_level', 'n_samples_variant_level', 'samplesize','n_obs',
'num_cases','num_controls','n_effective_samplesize'
]

p_value_list = [
'p_value','p','pval','pvalue','p-value','p.metal','neg_log_10_p_value'
'p', 'pvalue', 'p_value', 'pval', 'p_val', 'gc_pvalue','hetpval',
'p2_value','p.value','acace_pval','ace_pval','ala_pval','alb_pval','apoa1_pval','apob_apoa1_pval','apob_pval','bohbut_pval','cit_pval','crea_pval',
'dha_fa_pval','dha_pval','estc_pval','faw3_fa_pval','faw3_pval','faw6_fa_pval','faw6_pval','freec_pval','glc_pval','gln_pval',
'glol_pval','gly_pval','gp_pval','hdl2_c_pval','hdl3_c_pval','hdl_c_pval','hdl_d_pval','hdl_tg_pval','his_pval','idl_c_pval',
'idl_ce_pval','idl_fc_pval','idl_l_pval','idl_p_pval','idl_pl_pval','idl_tg_pval','ile_pval','l_hdl_c_pval','l_hdl_ce_pval','l_hdl_fc_pval',
'l_hdl_l_pval','l_hdl_p_pval','l_hdl_pl_pval','l_hdl_tg_pval','l_ldl_c_pval','l_ldl_ce_pval','l_ldl_fc_pval','l_ldl_l_pval','l_ldl_p_pval','l_ldl_pl_pval',
'l_ldl_tg_pval','l_vldl_c_pval','l_vldl_ce_pval','l_vldl_fc_pval','l_vldl_l_pval','l_vldl_p_pval','l_vldl_pl_pval','l_vldl_tg_pval','la_fa_pval','la_pval',
'lac_pval','ldl_c_pval','ldl_d_pval','ldl_tg_pval','leu_pval','m_hdl_c_pval','m_hdl_ce_pval','m_hdl_fc_pval','m_hdl_l_pval','m_hdl_p_pval',
'm_hdl_pl_pval','m_hdl_tg_pval','m_ldl_c_pval','m_ldl_ce_pval','m_ldl_fc_pval','m_ldl_l_pval','m_ldl_p_pval','m_ldl_pl_pval','m_ldl_tg_pval','m_vldl_c_pval',
'm_vldl_ce_pval','m_vldl_fc_pval','m_vldl_l_pval','m_vldl_p_pval','m_vldl_pl_pval','m_vldl_tg_pval','mufa_fa_pval','mufa_pval','pc_pval','phe_pval',
'pufa_fa_pval','pufa_pval','pyr_pval','remnant_c_pval','s_hdl_c_pval','s_hdl_ce_pval','s_hdl_fc_pval','s_hdl_l_pval','s_hdl_p_pval','s_hdl_pl_pval',
's_hdl_tg_pval','s_ldl_c_pval','s_ldl_ce_pval','s_ldl_fc_pval','s_ldl_l_pval','s_ldl_p_pval','s_ldl_pl_pval','s_ldl_tg_pval','s_vldl_c_pval','s_vldl_ce_pval',
's_vldl_fc_pval','s_vldl_l_pval','s_vldl_p_pval','s_vldl_pl_pval','s_vldl_tg_pval','serum_c_pval','serum_tg_pval','sfa_fa_pval','sfa_pval','sm_pval',
'tg_pg_pval','totcho_pval','totfa_pval','totpg_pval','tyr_pval','unsat_pval','val_pval','vldl_c_pval','vldl_d_pval',
'vldl_tg_pval','xl_hdl_c_pval','xl_hdl_ce_pval','xl_hdl_fc_pval','xl_hdl_l_pval','xl_hdl_p_pval','xl_hdl_pl_pval','xl_hdl_tg_pval','xl_vldl_c_pval','xl_vldl_ce_pval',
'xl_vldl_fc_pval','xl_vldl_l_pval','xl_vldl_p_pval','xl_vldl_pl_pval','xl_vldl_tg_pval','xs_vldl_c_pval','xs_vldl_ce_pval','xs_vldl_fc_pval','xs_vldl_l_pval','xs_vldl_p_pval',
'xs_vldl_pl_pval','xs_vldl_tg_pval','xxl_vldl_c_pval','xxl_vldl_ce_pval','xxl_vldl_fc_pval','xxl_vldl_l_pval','xxl_vldl_p_pval','xxl_vldl_pl_pval','xxl_vldl_tg_pval','gc.pvalue',
'p-val','gc.pvalue','pval_eas','pval_ind','pval_ira','jass_pval','univariate_min_pval','univariate_min_qval','mtag_pval','pval_dds',
'pval_fe','pval_q','pval_re','pval_re2','pval_uganda','pval_aadm','pval_dcc','pval_dds','pval_fe','pval_q',
'pval_re','pval_re2','pval_uganda','pval_dcc','pval_dds','pval_fe','pval_q','pval_re','pval_re2','pval_uganda',
'pvalue_all','het_p_value',
'european_ancestry_pval_fix','european_ancestry_pval_hetqtest','european_ancestry_pval_rand','multiancestry_pval_fix','multiancestry_pval_hetqtest','multiancestry_pval_rand',
'p-value_gc','_value','neg_log_10_p_value','meta.pval','het_pvalue','all.fixed.pval','all.random.pval','ff.fixed.pval','ff.random.pval','frequentist_add_pvalue',
'p-val','p-val','neg_log_10_p_value','score.pval','het_p_value','p-value_ancestry_het','p-value_association','p-value_residual_het','frequentist_add_pvalue','frequentist_add_pvalue',
'frequentist_add_pvalue','p-value_ancestry_het','p-value_association','p-value_residual_het','frequentist_add_wald_pvalue_1','p.value_ancestry_het','p.value_association','p.value_residual_het','robust_pva',
'robust_pval','meansphericaleqivalent:allmultiethnicmeta-analysisp-valuesforemmax-vt',
'p_score','freec_pval','p_bolt_lmm','p.value_association','p-value_association',
'pvalue_plink'
]

info_list =[
'info','info-score','info_all','info_score','info_ukbb','n_informative',
'impute_info','info.ukb','info_ukb','all_info','variant_info_score','additional_info',
'info.new','median_info','imputation','imputation_quality','imputationaccuracy','imputationinfo',
'impu_rsq','oncoarray_imputation_r2','imp_quality','imp_rsqr'
]

maf_list= ['frq1',
'eaf', 'frq','maf','frq_u', 'f_u',
'effect_allele_frequency','effectallelefreq','freq','freq1','a1freq',
'effect_allele_frequency_cases','effect_allele_frequency_controls',
'freq_effect_allele','eff_allele_freq','a1_freq','minorallelefreq.cases',
'minorallelefreq.controls','freq_a1','freqa1','ref_allele_frequency',
'freq.allele1.hapmapceu','allelefreq.cases','allelefreq.controls',
'minorallelefrequency','freq.allele1.hapmapceu','freq.a1.esp.eur','freq1.hapmap',
'freqa1','eff_all_freq','1000g_allele_freq','effect_allele_freq','freq_effect_allele_all',
'effectallelefrequencyeaf','effect_allele_freq','coded_allele_freq','effect_allle_frequency',
'effect_allele_fre','effect_allele_frequenc','effect_allele_freq','effect_allele_freque',
'_allele_frequency','effect_allele_freque','frequency','t_allele_frequen',
'freq_effect_allele','allele_frequencies','effect_allele_freq','all.freq.var',
'ff.freq.var','eff_allele_freq','effect-allele-frequency','freqa1','effect_allele_freq',
'allelefreq','eafreq','maxeafreq','mineafreq','expected_minor_allele_freq',
'freq.a1.1000g.eur','freq_hapmap','a1.af','af','af1',
'cases_maf', 'controls_maf', 'global_maf', 'cases_maf', 'controls_maf',
'global_maf', 'eaf_ukb', 'maf_ukb', 'eaf_a1', 'af_coded_all',
'eaf_eas', 'eaf_eur', 'eaf_ind', 'eaf_ira', 'maf_nw', 'af_dds',
'af_uganda', 'af_aadm', 'af_dcc', 'af_dds', 'af_uganda',
'af_dcc', 'af_dds', 'af_uganda', 'eaf_cases', 'eaf_controls',
'bcac_gwas_all_eaf_controls', 'bcac_gwas_erneg_eaf_controls',
'bcac_gwas_erpos_eaf_controls', 'bcac_icogs2_eaf_controls',
'bcac_icogs2_erneg_eaf_controls', 'bcac_icogs2_erpos_eaf_controls',
'bcac_onco2_eaf_controls', 'bcac_onco2_erneg_eaf_controls',
'bcac_onco2_erpos_eaf_controls', 'bcac_onco_icogs_gwas_eaf_controls',
'bcac_onco_icogs_gwas_erneg_eaf_controls',
'bcac_onco_icogs_gwas_erpos_eaf_controls',
'bcac_icogs1_european_controls_eaf', 'bcac_icogs_survival_eaf',
'bcac_icogs_survival_erneg_eaf', 'bcac_icogs_survival_erpos_eaf',
'eaf_hapmap_yri', 'eaf_hapmap_gih', 'eaf_hapmap_ceu',
'eaf_hapmap_chb', 'maf_hapmap_rel27', 'maf_co', 'cases_maf',
'controls_maf', 'avreaf', 'effect_af', 'aaf', 'af_coded_all',
'maf_ukb', 'all_meta_af', 'var1_maf', 'all_maf', 'eaf.cad',
'eaf.lipid', 'all_maf', 'all_maf', 'cases_maf', 'controls_maf',
'all_maf', 'eaf_dnk', 'eaf_fra', 'eaf_roi', 'exac_sas_maf',
'sas_maf', 'amr_maf', 'exac_amr_maf', 'af1',


]

beta_list = ['p_bolt_lmm_inf',
    'beta', 'b','effect','effects','n_effective','effect_a2',
    'effects','effect1','meta.effect','effect.1','effec','effec.1',
    'effect_a1','effect_all','effect_nw','a1_effect','maineffects',
    'frequentist_add_beta_1:add/plink_pheno=1','beta_plink','effect_size',
    'betazscale'
]
 
se_list = [
    'se', 'se_plink','se_saige','se_lm','standard_error',
    'standarderror','standard_error_of_beta','standard',
    'standard_error_asian','standard_error_black','standard_error_chinese', 
    'standard_error_white', 'standarderror', 'standard_error_site1',
    'standard_error_site2','se_error','tandard_error','dard_error',
    'andard_error','dard_error','tandard_error','_error','serror',
    'stderr','stderr_all','stderrlogor','stderr_all','stderr_nw',
    'stderr_females','stderr_males','log_or_ste','logor_se','or_se','log_odds_se',
    'logor.se','stderrlogor','log_or_ste','se1','se2','sebeta','mtag_se','freqse',
    'sebeta_snp_add','robust.se','se.coef.','se_gc','meta.se','sebeta','freqse',
    'odds_ratio_se','se_of_beta','se_0','se_1','se_2','se_3','freqse','frequentist_add_se_1',
    'freqse'
    ]

or_list =[
'or','log_odds','log_or','logor','or_random','odds_ratio','oddsratiominorallele',
'hm_odds_ratio','odds_rat','odd','odds_ration','odds','odds_ra','odds_rat','__odds_ratio__',
'heterozygous_or','homozygous_or'
]
 
zscore = [
  'z.weightedsumz', 'z_value','mtag_z','z.meta', 
  'gc.zscore','z_stat','zscore','z-score','z_score',
  'gz_zscore','gz-score','z','__z__','zval',
  'weighted_z','zscore_metal','betazscale'
]
direction = [
    'direction','direction_ukbb_tagc','direction_effects_cohorts','direction_effects_cohorts_all','direction_effects_cohorts_females','direction_effects_cohorts_males',
    'direct','direction_by_study','effect_direction','mantra_dir',
    'dir','cohort_dir' 
]

def count_unique_values(lst):
    # Split the list by comma and remove leading/trailing spaces
    values = lst
    # Find unique values
    unique_values = set(values)
    # Return the count of unique values
    return len(unique_values)

# Print the count of unique values for each list
print(f"Chromosome List: {count_unique_values(chromosome_list)} unique values")
print(f"Snp List: {count_unique_values(snp_list)} unique values")
print(f"Base Pair List: {count_unique_values(base_pair_list)} unique values")
print(f"Effect Allele List: {count_unique_values(effect_allele_list)} unique values")
print(f"Alternative Allele List: {count_unique_values(alternative_allele_list)} unique values")
print(f"N List: {count_unique_values(N_list)} unique values")
print(f"P Value List: {count_unique_values(p_value_list)} unique values")
print(f"Info List: {count_unique_values(info_list)} unique values")
print(f"Maf List: {count_unique_values(maf_list)} unique values")
print(f"Beta List: {count_unique_values(beta_list)} unique values")
print(f"Se List: {count_unique_values(se_list)} unique values")
print(f"Or List: {count_unique_values(or_list)} unique values")
print(f"Zscore: {count_unique_values(zscore)} unique values")
print(f"Direction: {count_unique_values(direction)} unique values")

exit(0)


def create_directory(directory):
    """Function to create a directory if it doesn't exist."""
    if not os.path.exists(directory):  # Checking if the directory doesn't exist
        os.makedirs(directory)  # Creating the directory if it doesn't exist
    return directory  # Returning the created or existing directory


def get_file_extension(file_path):
    _, file_extension = os.path.splitext(file_path)
    return file_extension


def detect_delimiter(file_path):
    import csv
    # Use csv.Sniffer to automatically detect the delimiter
    try:
        with open(file_path, 'r') as f:
            sample_data = f.read(1024)  # Read a sample of the file
            dialect = csv.Sniffer().sniff(sample_data)
        
        # Access the detected delimiter from the dialect
        detected_delimiter = dialect.delimiter
        return detected_delimiter
    except:
        return '\t'
 
    # Return the detected delimiter
    
def detect_quotes(file_path, lines_to_check=5):
    with open(file_path, 'r', encoding='utf-8') as f:
        # Read the first few lines
        sample_lines = [f.readline().strip() for _ in range(lines_to_check)]

    # Check if double quotes are consistently used as a separator
    double_quotes_used = any('"' in line for line in sample_lines)
    
    return double_quotes_used

def remove_quotes(gwaspath):
    with open(gwaspath, 'r', encoding='utf-8') as input_file:
        content = input_file.read()            
    modified_content = content.replace('"', '').replace("\t",",")
    
    with open(gwaspath+".modified", 'w', encoding='utf-8') as output_file:
        output_file.write(modified_content)
    return gwaspath+".modified"





def detect_delimitor_gz(gwaspath):
    from detect_delimiter import detect
    delimiter1 = ""
    with open(gwaspath) as f:
        content = '\n'.join(f.readline() for _ in range(20))
        delimiter1 = detect(content)
    delimiter2 = detect_delimiter(gwaspath)
 
    possible_delimiters = ['\t', ';', ',', '|','\s+',' ']
 
    def read_and_count_columns(file_path, delimiter):
        try:
            x = pd.read_csv(file_path,sep=delimiter, comment='#',encoding="utf-8")
            return x.shape[1]
        except:
            return 1
 
    bestdelimtor = ""
    for delimiter in possible_delimiters:
        columns = read_and_count_columns(gwaspath, delimiter)
 
        if columns < 2:
            continue
        else:
            bestdelimtor = delimiter
            
    
    if len(bestdelimtor) == 0:
        print("Cannot find delimiter",gwaspath)
        print(detect_quotes(gwaspath, lines_to_check=5))
        x = pd.read_csv(gwaspath, comment='#',encoding="utf-8")
 
        return ""
    else:
        return bestdelimtor
      
 
      
def list_and_sort_files_by_size(directory_path):
    try:
        # List files in the directory
        files = os.listdir(directory_path)

        # Sort the list based on file size
        sorted_files = sorted(files, key=lambda file: os.path.getsize(os.path.join(directory_path, file)))

        # Return the sorted filenames
        return sorted_files

    except FileNotFoundError:
        #print(f"Directory not found: {directory_path}")
        return None
 

def list_files_with_sizes(directory):
    file_info_list = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            file_size = os.path.getsize(file_path)
            file_info_list.append((file_path, file_size))
    return file_info_list

def find_largest_file(file_info_list):
    if not file_info_list:
        return None
    
    largest_file = max(file_info_list, key=lambda x: x[1])
    return largest_file



def rename_file_ending_with_digit(directory_path):
    if not os.path.exists(directory_path):
        print(f"Directory '{directory_path}' does not exist.")
        return

    # Iterate through all files in the directory
    for filename in os.listdir(directory_path):
        file_path = os.path.join(directory_path, filename)

        
        # Check if the file ends with a number
        if filename[-1].isdigit() and filename[-2]=="." and "gz" in filename:
            # Rename the file to "gz"
            new_file_path = os.path.join(directory_path, filename[:-2])
            os.rename(file_path, new_file_path)
            print(f"File renamed to: {new_file_path}")
        else:
            #print(f"Skipped file '{filename}' - does not end with a number.")
            pass


def replace_spaces_with_underscores(directory):
    # Iterate over all files and directories in the given directory
    for root, dirs, files in os.walk(directory):
        for file in files:
            old_path = os.path.join(root, file)
            new_path = os.path.join(root, file.replace(' ', '_'))
            os.rename(old_path, new_path)
        for dir_name in dirs:
            old_path = os.path.join(root, dir_name)
            new_path = os.path.join(root, dir_name.replace(' ', '_'))
            os.rename(old_path, new_path)
def poker(workingdirec,directory,url):
    sorted_files = list_and_sort_files_by_size(workingdirec+os.sep+directory)
    
    if sorted_files:
        filename = sorted_files[-1]
    else:
        return
     
    extension = get_file_extension(filename)

    if extension in ['.xlsx']:
        gwaspath = workingdirec+os.sep+directory+os.sep+filename
        create_directory(workingdirec+os.sep+directory+os.sep+"temp")
        
        try:
            gwasfile = pd.read_excel(gwaspath)
            gwaspath = workingdirec+os.sep+directory+os.sep+"gwas.csv"
            gwasfile.to_csv(gwaspath, index=False)
        except:
            print("Cannot read xlxs file! Kindly download the full file!")
            return

    if extension in [".zip"]:
        create_directory(workingdirec+os.sep+directory+os.sep+"temp")
        
        subprocess.run("bash -c '7z x "+workingdirec+os.sep+directory+os.sep+filename+" -o"+workingdirec+os.sep+directory+os.sep+"temp -y'", shell=True)
        replace_spaces_with_underscores(workingdirec+os.sep+directory+os.sep+"temp")

        tarextractedpath = workingdirec+os.sep+directory+os.sep+"temp"
        
        files_with_sizes = list_files_with_sizes(tarextractedpath)
        largest_file = find_largest_file(files_with_sizes)
        largest_file_path, largest_file_size = largest_file
        
        filename = largest_file_path.split("/")[3:]
        filename = "/".join(filename)
        extension = get_file_extension(largest_file_path)
        print(extension)
        print(filename)
        #exit(0)

    if extension in [".gz",".tar",".1"] and ".tar" in filename:
        create_directory(workingdirec+os.sep+directory+os.sep+"temp")
        
        if extension in [".tar"] and ".gz" not in filename:
            subprocess.run("bash -c 'tar -xvf "+workingdirec+os.sep+directory+os.sep+filename+" -C "+workingdirec+os.sep+directory+os.sep+"temp'", shell=True)
            replace_spaces_with_underscores(workingdirec+os.sep+directory+os.sep+"temp")

        else:
            subprocess.run("bash -c 'tar -xvzf "+workingdirec+os.sep+directory+os.sep+filename+" -C "+workingdirec+os.sep+directory+os.sep+"temp'", shell=True)
            replace_spaces_with_underscores(workingdirec+os.sep+directory+os.sep+"temp")

        tarextractedpath = workingdirec+os.sep+directory+os.sep+"temp"
        
        files_with_sizes = list_files_with_sizes(tarextractedpath)
 
        largest_file = find_largest_file(files_with_sizes)
        largest_file_path, largest_file_size = largest_file
        
        filename = largest_file_path.split("/")[3:]
        filename = "/".join(filename)
        extension = get_file_extension(largest_file_path)
        gwaspath = workingdirec+os.sep+directory+os.sep+filename
        try:
            gwasfile = pd.read_csv(gwaspath,nrows=20,compression='gzip', comment='#',encoding="utf-8")
        except:
            print("Cannot read gzip file as gz",gwaspath)
            print("Trying as normal csv")
            try:
                gwasfile = pd.read_csv(gwaspath,nrows=20, comment='#',encoding="utf-8")
                
            except Exception as e:
                if "EmptyDataError" in type(e).__name__ or "ParserError" in type(e).__name__:
                    print("EmptyDataFrame or ParserError!")
                    
                print("Trying with gunzip")
                
                try:
                    subprocess.run("bash -c 'gunzip -c "+workingdirec+os.sep+directory+os.sep+filename+" | head -n 100 > "+workingdirec+os.sep+directory+os.sep+"temp.csv'", shell=True)
                    #print("gunzip -c "+"allgwas"+os.sep+directory+os.sep+filename+" > "+"allgwas"+os.sep+directory+os.sep+"temp.csv")
                    #exit(0)
                    gwasfile = pd.read_csv(workingdirec+os.sep+directory+os.sep+"temp.csv",nrows=20, comment='#',encoding="utf-8")
                    
                except Exception as e:
                    if "ParserError" in type(e).__name__:
                        print("ParserError")
                        return
                    try:
                        subprocess.run("bash -c 'cat "+workingdirec+os.sep+directory+os.sep+filename+" > "+workingdirec+os.sep+directory+os.sep+"temp.csv'")
                        gwasfile = pd.read_csv(workingdirec+os.sep+directory+os.sep+"temp.csv",nrows=20, comment='#',encoding="utf-8")
                    
                    except:
                        print("Trying with double gunzip!")
                        try:
                            subprocess.run("bash -c 'gunzip -c "+workingdirec+os.sep+directory+os.sep+filename+" | gunzip -c | head -n 100 > "+workingdirec+os.sep+directory+os.sep+"temp.csv'", shell=True)
                            gwasfile = pd.read_csv(workingdirec+os.sep+directory+os.sep+"temp.csv",nrows=20, comment='#',encoding="utf-8")
                        
                        except:
                            print("Trying with zcat!")
                            
                            try:
                                subprocess.run("zcat  "+workingdirec+os.sep+directory+os.sep+filename+" | head -n 100 > "+workingdirec+os.sep+directory+os.sep+"temp.csv")
                                gwasfile = pd.read_csv(workingdirec+os.sep+directory+os.sep+"temp.csv",nrows=20, comment='#',encoding="utf-8")
                                print(gwasfile.head())
                            except:
                                        
                                print("CANNOT READ!")
                                exit(0)




    if extension in [".tsv",".txt",".csv",".1_summary_table",".ma",".assoc",".meta",".tbl",".linear",".logistic"]:
        gwaspath = workingdirec+os.sep+directory+os.sep+filename
        #print(workingdirec,directory,filename)
        #exit(0)
        try:
            gwasfile = pd.read_csv(gwaspath,nrows=20, comment='#',encoding="utf-8")
        
        except: 
            print("Trying with cat ")
            try:
                #print(workingdirec+os.sep+directory+os.sep+filename)
                subprocess.run("bash -c 'cat  "+workingdirec+os.sep+directory+os.sep+filename+" | head -n 100 > "+workingdirec+os.sep+directory+os.sep+"temp.csv'", shell=True)
                gwasfile = pd.read_csv(workingdirec+os.sep+directory+os.sep+"temp.csv",nrows=20, comment='#',encoding="utf-8")
                #print(gwasfile.head())
            except Exception as e:
                if "EmptyDataError" in type(e).__name__ or "ParserError" in type(e).__name__:
                    print("EmptyDataFrame or ParserError")
                    return
            gwasfile = pd.read_csv(gwaspath,nrows=20, comment='#',encoding="utf-8")
        gwasfile.to_csv(workingdirec+os.sep+directory+os.sep+"gwas.csv",index=False)

    
    if extension in [".gz",".GZ",".1"] and ".tar" not in filename:
        #print("allgwas"+os.sep+directory+os.sep+filename)
        gwaspath = workingdirec+os.sep+directory+os.sep+filename
        try:
            gwasfile = pd.read_csv(gwaspath,nrows=20,compression='gzip', comment='#',encoding="utf-8")
            #print(gwasfile.head())
        except:
            #print(";")
            try:
                gwasfile = pd.read_csv(gwaspath,nrows=20, comment='#',encoding="utf-8")
           
            except: 
                print("Trying with gunzip")
                
                try:
                    subprocess.run("bash -c 'gunzip -c "+workingdirec+os.sep+directory+os.sep+filename+" | head -n 100 > "+workingdirec+os.sep+directory+os.sep+"temp.csv'", shell=True)
                    gwasfile = pd.read_csv(workingdirec+os.sep+directory+os.sep+"temp.csv",nrows=20, comment='#',encoding="utf-8")
                except Exception as e:
                    if "ParserError" in type(e).__name__:
                        return

                    try:
                        subprocess.run("bash -c 'cat "+workingdirec+os.sep+directory+os.sep+filename+" > "+workingdirec+os.sep+directory+os.sep+"temp.csv'")
                        gwasfile = pd.read_csv(workingdirec+os.sep+directory+os.sep+"temp.csv",nrows=20, comment='#',encoding="utf-8")
                    except:
                        print("Trying with double gunzip!")
                        try:
                            subprocess.run("bash -c 'gunzip -c "+workingdirec+os.sep+directory+os.sep+filename+" | gunzip -c | head -n 100 > "+workingdirec+os.sep+directory+os.sep+"temp.csv'", shell=True)
                            gwasfile = pd.read_csv(workingdirec+os.sep+directory+os.sep+"temp.csv",nrows=20, comment='#',encoding="utf-8")
                            print(gwasfile.head())
                        except:
                            #subprocess.run("zcat  "+"allgwas"+os.sep+directory+os.sep+filename+" | head -n 100 > "+"allgwas"+os.sep+directory+os.sep+"temp.csv")
                            
                            try:
                                subprocess.run("zcat  "+workingdirec+os.sep+directory+os.sep+filename+" | head -n 100 > "+workingdirec+os.sep+directory+os.sep+"temp.csv")
                                gwasfile = pd.read_csv(workingdirec+os.sep+directory+os.sep+"temp.csv",nrows=20, comment='#',encoding="utf-8")
                            except:
                                        
                                print("CANNOT READ!")
                                return
        
        gwasfile.to_csv(workingdirec+os.sep+directory+os.sep+"gwas.csv",index=False)
        gwaspath = workingdirec+os.sep+directory+os.sep+"gwas.csv"

        newpath = ""
        if detect_quotes(gwaspath, lines_to_check=30):
            newpath = remove_quotes(gwaspath)
            print(newpath)
        else:
            #print("DEL")
            newpath = gwaspath
        bestdelimtor = detect_delimitor_gz(newpath)
        



        gwasfile = pd.read_csv(gwaspath, comment='#',sep=bestdelimtor,encoding="utf-8")
        gwasfile.to_csv(workingdirec+os.sep+directory+os.sep+"gwas.csv",index=False)
        
    return filename
        
    


def scrape(workingdirec,url,direc):
    def convert_to_numeric(size_str):
        if size_str=="-":
            return 0
        if size_str[-1] not in ['K','M','G','T']:
            return 0
        unit = size_str[-1]
        value = float(size_str[:-1])
        return value * size_mapping[unit]

    url = url
    response = requests.get(url)
    #print("Processing!",url)
    data = []
    
    if response.status_code == 200:
        soup = BeautifulSoup(response.text, 'html.parser')

        rows = soup.find_all('tr')[2:]
        rows = rows[:len(rows)-1]

        for row in rows[1:]:  # Skip the first row as it contains header information
            columns = row.find_all(['td', 'th'])
            name = columns[1].text.strip()
            last_modified = columns[2].text.strip()
            size = columns[3].text.strip()
            description = columns[4].text.strip()

            # Extracting URL if available
            if 'href' in columns[1].a.attrs:
                file_url = url +"/"+ columns[1].a['href']
            data.append([name, last_modified, size, description, file_url])
    else:
        print("URL does not exist",url)
        return
    
    columns = ["Name", "Last modified", "Size", "Description","URL"]

    df = pd.DataFrame(data, columns=columns)
    
    size_mapping = {'K': 1e3, 'M': 1e6, 'G': 1e9, 'T': 1e12}
    # Add a new column with numeric values
    df['Sizes_numeric'] = df['Size'].apply(convert_to_numeric)
    # Sort the DataFrame based on the 'Sizes_numeric' column
    df_sorted = df.sort_values(by='Sizes_numeric')
     
    max_index = df['Sizes_numeric'].idxmax()
    row_with_max_value = df.loc[max_index]

    finalurl = row_with_max_value["URL"]
    filename = finalurl.split("/")[-1]
    directoryname = finalurl.split("/")[-2]
    
    finalnameofthefile = row_with_max_value["Name"]
    
    if finalnameofthefile.endswith("/"):
        entries_ending_with_slash = df_sorted[df_sorted["Name"].str.endswith('/')]["Name"].values   
        if "harmonised/" in df_sorted['Name'].values:
            newurl = url = url+"/harmonised/"            
            response = requests.get(url)
            if response.status_code == 200:
                soup = BeautifulSoup(response.text, 'html.parser')

                rows = soup.find_all('tr')[2:]
                rows = rows[:len(rows)-1]

                for row in rows[1:]:  
                    columns = row.find_all(['td', 'th'])
                    name = columns[1].text.strip()
                    last_modified = columns[2].text.strip()
                    size = columns[3].text.strip()
                    description = columns[4].text.strip()

                    if 'href' in columns[1].a.attrs:
                        file_url = url +"/"+ columns[1].a['href']
                    
                    data.append([name, last_modified, size, description, file_url])
            columns = ["Name", "Last modified", "Size", "Description","URL"]

            df = pd.DataFrame(data, columns=columns)
            size_mapping = {'K': 1e3, 'M': 1e6, 'G': 1e9, 'T': 1e12}
            df['Sizes_numeric'] = df['Size'].apply(convert_to_numeric)

            # Sort the DataFrame based on the 'Sizes_numeric' column
            df_sorted = df.sort_values(by='Sizes_numeric')
                        
            max_index = df['Sizes_numeric'].idxmax()
            row_with_max_value = df.loc[max_index]

            finalurl = row_with_max_value["URL"]
            filename = finalurl.split("/")[-1]
            directoryname = finalurl.split("/")[-2]
            finalurl = finalurl.replace("%20"," ")
            filename = filename.replace("%20"," ")
            os.system("timeout -s KILL 10 wget -q "+"\"" +finalurl+"\""+ " -P "+workingdirec+os.sep+direc+os.sep)
            
            os.system("rm -rf wget-log*")
    
            filtered_rows = df[df['Name'].str.contains("readme", case=False)]
            if len(filtered_rows)>0:
                #print(filtered_rows)
                uu = filtered_rows["URL"].values[0]
                os.system("wget -q "+uu+" -P "+workingdirec+os.sep+direc+os.sep)
    
            df_sorted = df_sorted.drop('Sizes_numeric', axis=1)
            df_sorted.to_csv(workingdirec+os.sep+direc+os.sep+"Information.csv")
        else:
            entries_ending_with_slash = df_sorted[df_sorted["Name"].str.endswith('/')]["Name"].values   
            for loop in entries_ending_with_slash:
                newurl = url = url+os.sep+loop           
                response = requests.get(url)
                if response.status_code == 200:
                    soup = BeautifulSoup(response.text, 'html.parser')
                    rows = soup.find_all('tr')[2:]
                    rows = rows[:len(rows)-1]
                    for row in rows[1:]:  # Skip the first row as it contains header information
                        columns = row.find_all(['td', 'th'])
                        name = columns[1].text.strip()
                        last_modified = columns[2].text.strip()
                        size = columns[3].text.strip()
                        description = columns[4].text.strip()

                        if 'href' in columns[1].a.attrs:
                            file_url = url +"/"+ columns[1].a['href']
                        
                        data.append([name, last_modified, size, description, file_url])

                columns = ["Name", "Last modified", "Size", "Description","URL"]

                df = pd.DataFrame(data, columns=columns)
                size_mapping = {'K': 1e3, 'M': 1e6, 'G': 1e9, 'T': 1e12}
                def convert_to_numeric(size_str):
                    if size_str=="-":
                        return 0
                    if size_str[-1] not in ['K','M','G','T']:
                        return 0
                    unit = size_str[-1]
                    value = float(size_str[:-1])
                    return value * size_mapping[unit]

                df['Sizes_numeric'] = df['Size'].apply(convert_to_numeric)

                df_sorted = df.sort_values(by='Sizes_numeric')
                max_index = df['Sizes_numeric'].idxmax()
                row_with_max_value = df.loc[max_index]

                finalurl = row_with_max_value["URL"]
                filename = finalurl.split("/")[-1]
                directoryname = finalurl.split("/")[-2]
                
                finalurl = finalurl.replace("%20"," ")
                filename = filename.replace("%20"," ")
                os.system("timeout -s KILL 10 wget -q "+"\"" +finalurl+"\""+ " -P "+workingdirec+os.sep+direc+os.sep)
                
                os.system("rm -rf wget-log*")
        
                filtered_rows = df[df['Name'].str.contains("readme", case=False)]
                if len(filtered_rows)>0:
                    #print(filtered_rows)
                    uu = filtered_rows["URL"].values[0]
                    os.system("wget -q "+uu+" -P "+workingdirec+os.sep+direc+os.sep)
        
                df_sorted = df_sorted.drop('Sizes_numeric', axis=1)
                df_sorted.to_csv(workingdirec+os.sep+direc+os.sep+"Information.csv")    
    else:
        finalurl = finalurl.replace("%20"," ")
        filename = filename.replace("%20"," ")
        
        os.system("timeout -s KILL 10 wget -q "+"\"" +finalurl+"\""+ " -P "+workingdirec+os.sep+direc+os.sep)
        os.system("rm -rf wget-log*")
 
        filtered_rows = df[df['Name'].str.contains("readme", case=False)]
        if len(filtered_rows)>0:
            #print(filtered_rows)
            uu = filtered_rows["URL"].values[0]
            os.system("wget -q "+uu+" -P "+workingdirec+os.sep+direc+os.sep)
 
        df_sorted = df_sorted.drop('Sizes_numeric', axis=1)
        df_sorted.to_csv(workingdirec+os.sep+direc+os.sep+"Information.csv")
 
 
from IPython.core.display import HTML

# Include the Bootstrap CSS link
bootstrap_css_link = '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css" integrity="sha384-Gn5384xqQ1aoWXA+058RXPxPg6fy4IWvTNh0E263XmFcJlSAwiGgFAW/dAiS6JXm" crossorigin="anonymous">'
HTML(bootstrap_css_link)

def remove_quotes2(gwaspath):
    with open(gwaspath, 'r', encoding='utf-8') as input_file:
        content = input_file.read()            
    modified_content = content.replace('"',"").replace(":","_").replace("\t",",")
    
    with open(gwaspath+".modified", 'w', encoding='utf-8') as output_file:
        output_file.write(modified_content)
    return gwaspath+".modified"

def detect_delimitor_gz2(gwaspath):
    from detect_delimiter import detect
    delimiter1 = ""
    with open(gwaspath) as f:
        content = '\n'.join(f.readline() for _ in range(20))
        delimiter1 = detect(content)
    delimiter2 = detect_delimiter(gwaspath)
 
    possible_delimiters = [',','\s+']
 
    def read_and_count_columns(file_path, delimiter):
        try:
            x = pd.read_csv(file_path,sep=delimiter, comment='#',encoding="utf-8")
            return x.shape[1]
        except:
            return 1

        

    # Iterate over each delimiter
    bestdelimtor = ""
    for delimiter in possible_delimiters:
        columns = read_and_count_columns(gwaspath, delimiter)
        #print(columns)
        if columns < 2:
            continue
        else:
            bestdelimtor = delimiter
            
    
    if len(bestdelimtor) == 0:
        #print("Cannot find delimiter",gwaspath)
        #print(detect_quotes(gwaspath, lines_to_check=5))
        #x = pd.read_csv(gwaspath, comment='#',encoding="utf-8")
        #print(x.head())
        #exit(0)
        return "\s+"
    else:
        return bestdelimtor

from collections import Counter

def is_string_in_list(string_list, target_list):
    finallist = []
    for item in string_list:
        if item in target_list:
            finallist.append(item)

    return finallist
def doi_to_bibtex(doi):
    base_url = f'https://doi.org/{doi}'
    headers = {'Accept': 'application/x-bibtex'}

    try:
        response = requests.get(base_url, headers=headers)
        if response.status_code == 200:
            return response.text
        else:
            print(f"Error: {response.status_code}")
    except Exception as e:
        print("An error occurred:", e)
    
    return None

def pmid_to_doi(pmid):
    base_url = f'https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids={pmid}&format=json'
    try:
        response = requests.get(base_url)
        data = response.json()
        records = data.get('records', [])
        if records:
            for record in records:
                doi = record.get('doi')
                if doi:
                    return doi
    except Exception as e:
        print("An error occurred:", e)

    return None
if __name__ == "__main__":
    # Create an ArgumentParser
    #print("Download wget for windows and place it in the same directory: ","https://eternallybored.org/misc/wget/")
    
    parser = argparse.ArgumentParser(description='List GWAS Columns.')
    parser.add_argument('--processedfile', type=str, help='Specify the GWAS file path.')
    
    args = parser.parse_args()
 
    gwas_file = args.processedfile
    

    create_directory(gwas_file.split(".")[0].strip().replace(" ","_"))
    workingdirec = gwas_file.split(".")[0].strip().replace(" ","_")

    workingdirec = workingdirec+os.sep+"allgwas"
    create_directory(workingdirec)
    print(gwas_file)
    data = pd.read_csv(gwas_file,encoding='utf-8') 
    
    print("Searching for ",len(pd.read_csv(gwas_file))," on GWAS catalog!")
    print("Data will be stored in :",workingdirec)
    
    # Call the function to process the GWAS file
    df = pd.read_csv(gwas_file,index_col=0)
    print("\n")
    print("Processing the following file!")
    print(df.head())
    print("\n")
    
    
    import sys
    df = df.reset_index()
    urls = df['summaryStatistics'].values
    directories = df['accessionId'].values
    traits = df['reportedTrait'].values
    pubmedid = df["pubmedId"].values
    #print(pubmedid[0])
    import requests
    

    # Example usage:


    for loop in range(0,len(urls)):
        create_directory(workingdirec+os.sep+directories[loop])
        files_in_directory = os.listdir(workingdirec+os.sep+directories[loop])
        print(f"Processing {traits[loop]}: downloading information from {urls[loop]} and saving it to {directories[loop]}")
                
        try:
            # Check if Information file exist!
            files_in_directory.remove("Information.csv")
        except:
            scrape(workingdirec,urls[loop],directories[loop])
            files_in_directory = os.listdir(workingdirec+os.sep+directories[loop])
            try:
                files_in_directory.remove("Information.csv")
            except:
                pass

        number_of_files = len(files_in_directory)
        
        if number_of_files>1:
            pass
        else:
            scrape(workingdirec,urls[loop],directories[loop])
            pass    
        
    for loop in range(0,len(directories)):
        rename_file_ending_with_digit(workingdirec+os.sep+directories[loop])  

    actualfilenames = []
    for loop in range(0,len(directories)):
        print(f"Normalizing extensions! {directories[loop]}: {urls[loop]}")
        actualfilenames.append(poker(workingdirec,directories[loop],urls[loop]))
    #exit(0)
    
    bootstrap_css_link = '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css" integrity="sha384-Gn5384xqQ1aoWXA+058RXPxPg6fy4IWvTNh0E263XmFcJlSAwiGgFAW/dAiS6JXm" crossorigin="anonymous">'

    # The rest of your code remains the same with slight modifications

    heading1 = ""
    for loop in range(0, len(directories)):

        directory_path = workingdirec + os.sep + directories[loop]

        heading1 += f"<div class='container mt-4 bg-light p-4 text-center'>"  # Added 'text-center' class here
        heading1 += f"<h1 class='text-primary'>Phenotype: {traits[loop]}</h1>"
        heading1 += f"<p class='text-muted'><strong>URL:</strong> {urls[loop]}</p>"


        try:
            gwaspath = directory_path + os.sep + "gwas.csv"
            remove_quotes2(gwaspath)
        except:
            pass

        try:
            gwaspath = directory_path + os.sep + "gwas.csv.modified"
            bestdelimtor = detect_delimitor_gz2(gwaspath)
            gwasfile = pd.read_csv(gwaspath, sep=",", comment='#', encoding="utf-8")
            if len(gwasfile.columns) < 2:
                gwasfile = pd.read_csv(gwaspath, sep=" ", comment='#', encoding="utf-8")

        except:
            gwasfile = "not"
            # print("File not found!")
            # continue
        row_data = pd.DataFrame(df.loc[loop]).transpose()
        heading1 += row_data.to_html(index=False, classes='table table-bordered')

        try:
            heading1 += "<h2 class='text-muted'>Information file</h2>"
            heading1 += pd.read_csv(workingdirec + os.sep + directories[loop] + os.sep + "Information.csv").to_html(
                index=False, classes='table table-bordered')
        except:
            heading1 += "<h2 class='text-muted'>Information file</h2>"
            heading1 += "<p>Information file is missing!</p>"

        if isinstance(gwasfile, str):
            heading1 += "<h2 class='text-muted'>GWAS file</h2>"
            heading1 += "<p>We were unable to process the GWAS file!</p>"
            continue

        gwasfile = gwasfile.head(2)
        gwasfile.columns = gwasfile.columns.str.lower().tolist()

        heading1 += "<h2 class='text-muted'>Processed file</h2>"
        heading1 += "<p>Processing file: " + actualfilenames[loop] + "!</p>"

        try:
            heading1 += "<h2 class='text-muted'>GWAS file head</h2>"
            heading1 += gwasfile.head().to_html(index=False, classes='table table-bordered')
        except:
            heading1 += "<h2 class='text-muted'>GWAS file head</h2>"
            heading1 += "<p>GWAS file was not downloaded, kindly do it manually!</p>"

        heading1 += "</div>"  # Closing div for styling

        col1 = is_string_in_list(gwasfile.columns.str.lower().tolist(), chromosome_list)
        col2 = is_string_in_list(gwasfile.columns.str.lower().tolist(), snp_list)
        col3 = is_string_in_list(gwasfile.columns.str.lower().tolist(), base_pair_list)
        col4 = is_string_in_list(gwasfile.columns.str.lower().tolist(), effect_allele_list)
        col5 = is_string_in_list(gwasfile.columns.str.lower().tolist(), alternative_allele_list)
        col6 = is_string_in_list(gwasfile.columns.str.lower().tolist(), N_list)
        col8 = is_string_in_list(gwasfile.columns.str.lower().tolist(), p_value_list)
        col9 = is_string_in_list(gwasfile.columns.str.lower().tolist(), info_list)
        col10 = is_string_in_list(gwasfile.columns.str.lower().tolist(), maf_list)
        col11 = is_string_in_list(gwasfile.columns.str.lower().tolist(), beta_list)
        col12 = is_string_in_list(gwasfile.columns.str.lower().tolist(), se_list)
        col13 = is_string_in_list(gwasfile.columns.str.lower().tolist(), or_list)
        col14 = is_string_in_list(gwasfile.columns.str.lower().tolist(), zscore)
        col15 = is_string_in_list(gwasfile.columns.str.lower().tolist(), direction)

        variable_list_mapping = {
            'col1': 'CHR',
            'col2': 'SNP',
            'col3': 'BP',
            'col4': 'A1',
            'col5': 'A2',
            'col6': 'N',
            'col8': 'P',
            'col9': 'INFO',
            'col10': 'MAF',
            'col11': 'BETA',
            'col12': 'SE',
            'col13': 'OR',
            'col14': 'Z',
            'col15': 'D'
        }
        allcolumns = []
        heading1 += f"<div class='container mt-4 bg-light p-4 text-center'>"  # Added 'text-center' class here

        # Your existing code
        #heading1 += "<h2 class='text-muted'>Original GWAS file headers</h2>"
        #heading1 += gwasfile.to_html(index=False, classes='table table-bordered')

        for variable, lst in variable_list_mapping.items():
            values = globals().get(variable)
            allcolumns.extend(values)
            print(f"{lst} -> {values}")
            heading1 += f"<p style='font-weight: bold; color: blue;'>{lst} -> {values}</p>"

        gwascols = gwasfile.columns.str.lower().tolist()
        unidentifiedcols = set(gwascols) - set(allcolumns)

        heading1 += "<h2 class='text-muted'>GWAS columns</h2>"
        heading1 += ', '.join(gwascols)
        heading1 += "<h2 class='text-muted'>Identified columns</h2>"
        heading1 += ', '.join(allcolumns)
        heading1 += "<h2 class='text-muted'>Unidentified Columns headers</h2>"
        heading1 += ', '.join(unidentifiedcols)


        pmid = pubmedid[loop]      
        doi = pmid_to_doi(pubmedid[loop])
 
        bibtexs = ""

        if doi:
            bibtex = doi_to_bibtex(doi)
            if bibtex:
                bibtexs = bibtex
                #print("BibTeX entry:")
                #print(bibtex)
            else:
                bibtexs = "Could not find bibtex"
        else:
            doi = "Could not find DOI!"

        heading1 += "<h2 class='text-muted'>PMID</h2>"
        heading1 += str(pmid)

        heading1 += "<h2 class='text-muted'>DOI</h2>"
        heading1 += doi

        bibtex_style = "font-family: monospace; white-space: pre-wrap;"

        # Constructing the HTML
        heading1 += "<div class='bibtex-container'>"  # Opening div for BibTeX
        heading1 += "<h2 class='text-muted'>BibTeX</h2>"
        heading1 += f"<div style='{bibtex_style}'>{bibtexs}</div>"
        heading1 += "</div>"  # Closing div for BibTeX container

        heading1 += "<br><br><br><br><br><br>"  # Additional spacing for better readability

        heading1 += "</div>"  # Closing div for styling
   
    #exit(0)
    # Save the generated HTML content to a file
    with open( gwas_file.split(".")[0].strip().replace(" ","_")+"_output.html", "w") as html_file:
        html_file.write(bootstrap_css_link + heading1)

    # Open the file in a web browser
    HTML(filename=gwas_file.split(".")[0].strip().replace(" ","_")+"_output.html")
