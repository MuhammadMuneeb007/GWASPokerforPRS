import pandas as pd
import numpy as np
import sys
import os
import argparse
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

def create_directory(directory):
    """Function to create a directory if it doesn't exist."""
    if not os.path.exists(directory):  # Checking if the directory doesn't exist
        os.makedirs(directory)  # Creating the directory if it doesn't exist
    return directory  # Returning the created or existing directory


def detect_quotes(file_path, lines_to_check=5):
    with open(file_path, 'r', encoding='utf-8') as f:
        # Read the first few lines
        sample_lines = [f.readline().strip() for _ in range(lines_to_check)]

    # Check if double quotes are consistently used as a separator
    double_quotes_used = any('"' in line for line in sample_lines)
    
    return double_quotes_used



def detect_delimitor_gz(gwaspath):
    from detect_delimiter import detect
   
 
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

        return ""
    else:
        return bestdelimtor
      



def remove_quotes(gwaspath):
    with open(gwaspath, 'r', encoding='utf-8') as input_file:
        content = input_file.read()            
    modified_content = content.replace('"', '').replace("\t",",")
    
    with open(gwaspath+".modified", 'w', encoding='utf-8') as output_file:
        output_file.write(modified_content)
    return gwaspath+".modified"


def get_file_extension(file_path):
    _, file_extension = os.path.splitext(file_path)
    return file_extension

def poker(workingdirec,directory,url):
    filename = directory
    directory = ""
    extension = get_file_extension(filename)

    if extension in ['.xlsx']:
        gwaspath = workingdirec+os.sep+directory+filename
        try:
            gwasfile = pd.read_excel(gwaspath)
            gwaspath = workingdirec+os.sep+directory+"gwas.csv"
            gwasfile.to_csv(gwaspath, index=False)
        except:
            print("Cannot read xlxs file! Kindly check file ",gwaspath)
            return

    if extension in [".zip"]:
        create_directory(workingdirec+os.sep+directory+os.sep+"temp")
        
        subprocess.run("bash -c '7z x "+workingdirec+os.sep+directory+filename+" -o"+workingdirec+os.sep+directory+"temp -y'", shell=True)
        replace_spaces_with_underscores(workingdirec+os.sep+directory+"temp")

        tarextractedpath = workingdirec+os.sep+directory+"temp"
        
        files_with_sizes = list_files_with_sizes(tarextractedpath)
        largest_file = find_largest_file(files_with_sizes)
        largest_file_path, largest_file_size = largest_file
        
        filename = largest_file_path.split("/")[1:]
        filename = "/".join(filename)
        extension = get_file_extension(largest_file_path)
        print(extension)
        print(filename)
        #exit(0)

    if extension in [".gz",".tar",".1"] and ".tar" in filename:
        create_directory(workingdirec+os.sep+directory+"temp")
        
        if extension in [".tar"] and ".gz" not in filename:
            subprocess.run("bash -c 'tar -xvf "+workingdirec+os.sep+directory+filename+" -C "+workingdirec+os.sep+directory+"temp'", shell=True)
            replace_spaces_with_underscores(workingdirec+os.sep+directory+"temp")

        else:
            subprocess.run("bash -c 'tar -xvzf "+workingdirec+os.sep+directory+filename+" -C "+workingdirec+os.sep+directory+"temp'", shell=True)
            replace_spaces_with_underscores(workingdirec+os.sep+directory+"temp")

        tarextractedpath = workingdirec+os.sep+directory+"temp"
        
        files_with_sizes = list_files_with_sizes(tarextractedpath)
 
        largest_file = find_largest_file(files_with_sizes)
        largest_file_path, largest_file_size = largest_file
        
        filename = largest_file_path.split("/")[1:]
        filename = "/".join(filename)
        extension = get_file_extension(largest_file_path)
        gwaspath = workingdirec+os.sep+directory+filename

        try:
            gwasfile = pd.read_csv(gwaspath,compression='gzip', comment='#',encoding="utf-8")
        except:
            print("Cannot read gzip file as gz",gwaspath)
            print("Trying as normal csv")
            try:
                gwasfile = pd.read_csv(gwaspath, comment='#',encoding="utf-8")
                
            except Exception as e:
                if "EmptyDataError" in type(e).__name__ or "ParserError" in type(e).__name__:
                    print("EmptyDataFrame or ParserError!")
                    
                print("Trying with gunzip")
                
                try:
                    subprocess.run("bash -c 'gunzip -c "+workingdirec+os.sep+directory+filename+"  > "+workingdirec+os.sep+directory+"temp.csv'", shell=True)
                    #print("gunzip -c "+"allgwas"+os.sep+directory+os.sep+filename+" > "+"allgwas"+os.sep+directory+os.sep+"temp.csv")
                    #exit(0)
                    gwasfile = pd.read_csv(workingdirec+os.sep+directory+"temp.csv",comment='#',encoding="utf-8")
                    
                except Exception as e:
                    if "ParserError" in type(e).__name__:
                        print("ParserError")
                        return
                    try:
                        subprocess.run("bash -c 'cat "+workingdirec+os.sep+directory+filename+" > "+workingdirec+os.sep+directory+"temp.csv'")
                        gwasfile = pd.read_csv(workingdirec+os.sep+directory+"temp.csv", comment='#',encoding="utf-8")
                    
                    except:
                        print("Trying with double gunzip!")
                        try:
                            subprocess.run("bash -c 'gunzip -c "+workingdirec+os.sep+directory+filename+" | gunzip -c  > "+workingdirec+os.sep+directory+"temp.csv'", shell=True)
                            gwasfile = pd.read_csv(workingdirec+os.sep+directory+"temp.csv",nrows=20, comment='#',encoding="utf-8")
                        
                        except:
                            print("Trying with zcat!")
                            
                            try:
                                subprocess.run("zcat  "+workingdirec+os.sep+directory+filename+" |  > "+workingdirec+os.sep+directory+"temp.csv")
                                gwasfile = pd.read_csv(workingdirec+os.sep+directory+"temp.csv", comment='#',encoding="utf-8")
                                print(gwasfile.head())
                            except:
                                        
                                print("CANNOT READ!")
                                exit(0)




    if extension in [".tsv",".txt",".csv",".1_summary_table",".ma",".assoc",".meta",".tbl",".linear",".logistic"]:
        gwaspath = workingdirec+os.sep+directory+filename
        #print(workingdirec,directory,filename)
        #exit(0)
        try:
            gwasfile = pd.read_csv(gwaspath, comment='#',encoding="utf-8")
        
        except: 
            print("Trying with cat ")
            try:
                #print(workingdirec+os.sep+directory+os.sep+filename)
                subprocess.run("bash -c 'cat  "+workingdirec+os.sep+directory+filename+" | > "+workingdirec+os.sep+directory+"temp.csv'", shell=True)
                gwasfile = pd.read_csv(workingdirec+os.sep+directory+"temp.csv",comment='#',encoding="utf-8")
                #print(gwasfile.head())
            except Exception as e:
                if "EmptyDataError" in type(e).__name__ or "ParserError" in type(e).__name__:
                    print("EmptyDataFrame or ParserError")
                    return
            gwasfile = pd.read_csv(gwaspath, comment='#',encoding="utf-8")
        gwasfile.to_csv(workingdirec+os.sep+directory+"gwas.csv",index=False)

    
   
    if extension in [".gz",".GZ",".1"] and ".tar" not in filename:
        #print("allgwas"+os.sep+directory+os.sep+filename)
        gwaspath = workingdirec+os.sep+directory+filename
        try:
            gwasfile = pd.read_csv(gwaspath,compression='gzip', comment='#',encoding="utf-8")
            #print(gwasfile.head())
        except:
            #print(";")
            try:
                gwasfile = pd.read_csv(gwaspath, comment='#',encoding="utf-8")
           
            except: 
                print("Trying with gunzip")
                
                try:
                    subprocess.run("bash -c 'gunzip -c "+workingdirec+os.sep+directory+filename+"  > "+workingdirec+os.sep+directory+"temp.csv'", shell=True)
                    gwasfile = pd.read_csv(workingdirec+os.sep+directory+"temp.csv", comment='#',encoding="utf-8")
                except Exception as e:
                    if "ParserError" in type(e).__name__:
                        return

                    try:
                        subprocess.run("bash -c 'cat "+workingdirec+os.sep+directory+filename+" > "+workingdirec+os.sep+directory+"temp.csv'")
                        gwasfile = pd.read_csv(workingdirec+os.sep+directory+"temp.csv", comment='#',encoding="utf-8")
                    except:
                        print("Trying with double gunzip!")
                        try:
                            subprocess.run("bash -c 'gunzip -c "+workingdirec+os.sep+directory+filename+" | gunzip -c |  > "+workingdirec+os.sep+directory+"temp.csv'", shell=True)
                            gwasfile = pd.read_csv(workingdirec+os.sep+directory+"temp.csv", comment='#',encoding="utf-8")
                        
                        except:
                            #subprocess.run("zcat  "+"allgwas"+os.sep+directory+os.sep+filename+" | head -n 100 > "+"allgwas"+os.sep+directory+os.sep+"temp.csv")
                            
                            try:
                                subprocess.run("zcat  "+workingdirec+os.sep+directory+filename+" |  > "+workingdirec+os.sep+directory+"temp.csv")
                                gwasfile = pd.read_csv(workingdirec+os.sep+directory+"temp.csv", comment='#',encoding="utf-8")
                            except:
                                        
                                print("CANNOT READ!")
                                return
        
        gwasfile.to_csv(workingdirec+os.sep+directory+"gwas.csv",index=False)
        gwaspath = workingdirec+os.sep+directory+"gwas.csv"

        newpath = ""
        if detect_quotes(gwaspath, lines_to_check=30):
            newpath = remove_quotes(gwaspath)
            
        else:
            #print("DEL")
            newpath = gwaspath

        bestdelimtor = detect_delimitor_gz(newpath)
        gwasfile = pd.read_csv(gwaspath, comment='#',sep=bestdelimtor,encoding="utf-8")
        gwasfile.to_csv(workingdirec+os.sep+directory+"gwas.csv",index=False)
        
    return filename

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

def remove_quotes2(gwaspath):
    with open(gwaspath, 'r', encoding='utf-8') as input_file:
        content = input_file.read()            
    modified_content = content.replace('"',"").replace(":","_").replace("\t",",")
    
    with open(gwaspath+".modified", 'w', encoding='utf-8') as output_file:
        output_file.write(modified_content)
    return gwaspath+".modified"


def detect_delimitor_gz2(gwaspath):
    from detect_delimiter import detect
  
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

        return "\s+"
    else:
        return bestdelimtor

def is_string_in_list(string_list, target_list):
    finallist = []
    for item in string_list:
        if item in target_list:
            finallist.append(item)

    return finallist
 
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='This module list the PRS columns in the GWAS file!')

    parser.add_argument('--gwasfile', type=str, help='Pass the GWAS file!')
    

    args = parser.parse_args()
    gwas_file = args.gwasfile
  

    data = pd.read_csv(gwas_file,encoding='utf-8') 
    print(data.head())
    
    import sys
    
    try:
        gwaspath = gwas_file
        remove_quotes2(gwaspath)
    except:
        pass
    
    

    try:
        gwaspath =  "gwas.csv.modified"
        bestdelimtor = detect_delimitor_gz2(gwaspath)
        gwasfile = pd.read_csv(gwaspath, sep=",", comment='#', encoding="utf-8")
        if len(gwasfile.columns) < 2:
            gwasfile = pd.read_csv(gwaspath, sep=" ", comment='#', encoding="utf-8")

    except:
        gwasfile = "not"
    
    
    if isinstance(gwasfile, str):
        print("Unable to extract GWAS. Kindly do it manually!")
        exit(0)
    
    
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

    transformeddic = {}
    originaldic = {}

    for variable, lst in variable_list_mapping.items():
        values = globals().get(variable)
        allcolumns.extend(values)
        print(f"{lst} -> {values}")
        originaldic[lst] = values

        if values:
            for v in values:
                transformeddic[v] = lst
        else:
            transformeddic["NA_"+lst] = lst

    print(transformeddic)

    for key, value in originaldic.items():
        print(key,"->",value)
    for key, value in transformeddic.items():
        print(key,"->",value)
            
    with open('transform.txt', 'w') as file:
        for key, value in transformeddic.items():
            file.write(f" {key} -> {value}\n")
            
    with open('transform1.txt', 'w') as file:
        for key, value in transformeddic.items():
            if "NA" in key:
                continue
            else:
                file.write(f" {key} -> {value}\n")

    from hugchat import hugchat
    from hugchat.login import Login
    sign = Login("YOUR HUGGING CHAT EMAIL", "YOUR HUGGING CHAT PASSWORD")
    cookies = sign.login()
    cookie_path_dir = "./cookies_snapshot"
    sign.saveCookiesToDir(cookie_path_dir)
    chatbot = hugchat.ChatBot(cookies=cookies.get_dict())  # or cookie_path="usercookies/<email>.json"
    
    file = 'transform.txt'
    import pandas as pd
    file_path = file
    try:
        with open(file_path, 'r') as file:
            content = file.read()
            print(content)
    except FileNotFoundError:
        print("File not found.")
    except IOError:
        print("Error reading the file.")


    file_path = 'Output.py'

    file = open(file_path, 'w')
    file.write("# ORGINAL MAPPING!\n\n")
    file.close()

    with open(file_path, 'a') as file:
        for key, value in originaldic.items():
            file.write(f"# {key} -> {value}\n")

    file = open(file_path, 'a')
    file.write("# Transformation for code !\n\n")
    file.close()
    
    with open(file_path, 'a') as file:
        for key, value in transformeddic.items():
            file.write(f"# {key} -> {value}\n")
            
    try:
        with open("transform1.txt", 'r') as file1:
            content = file1.read()
            
    except FileNotFoundError:
        print("File not found.")
    except IOError:
        print("Error reading the file.")

    query_result = chatbot.query("Convert the following transformation to pandas in python vertically and just give me mapping vertically for gwas.csv.modified: and output is finalgwas.csv "+"\n\n"+content)
    print(query_result)

    

    with open(file_path, 'a',encoding="utf-8") as file:
        file.write(query_result["text"])







        



