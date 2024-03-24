import pandas as pd
from fuzzywuzzy import fuzz
import requests
from bs4 import BeautifulSoup
import pandas as pd
import numpy as np
import pandas as pd
import numpy as np

import os
from collections import Counter
import pandas as pd
import numpy as np
import json
#from scholarly import scholarly


import os
import sys


from sklearn.feature_extraction.text import TfidfVectorizer
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import numpy as np
import pandas as pd
import nltk
import re
import os
import codecs
from sklearn import feature_extraction
#import mpld3
import pandas as pd
import re
from sklearn.decomposition import PCA
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import sys
from transformers import pipeline
from transformers import AutoTokenizer, AutoModelForQuestionAnswering
 
import io
def extract_number(input_string):
    try:
        # Use regular expression to find the first number in the string
        input_string = input_string.replace(",","").replace(".","")
        numbers = re.findall(r'\d+', input_string)
        
        # If numbers are found, return the first one as an integer; otherwise, return None
        return int(numbers[0].replace(',', '')) if numbers else None
    except:
        return 0

def questionanswer(q,textdata):
    try:    
        model_name='ahotrod/electra_large_discriminator_squad2_512'
        qa_pipeline = pipeline('question-answering', model=model_name, tokenizer=model_name)
        result = qa_pipeline(question=q, context=textdata)
        return result 
    except:
        return "-"
    
def getmeanswers_cases_controls(row,question):
    try:
        if "case" in str(row).lower() and "control" in str(row).lower():
            #print(question,row)
            #print(questionanswer(question,row)["answer"])
            #print(extract_number(questionanswer(question,row)["answer"]))
            return extract_number(questionanswer(question,row)["answer"])
        else:
            return "-"
    except:
        return "-"
    
def getmeanswers_samples(row,question):
    try:
        if "case" in str(row).lower() and "control" in str(row).lower():
            return "-"  
        else:
            #print(question,row)
            #print(questionanswer(question,row)["answer"])
            #print(extract_number(questionanswer(question,row)["answer"]))
            return extract_number(questionanswer(question,row)["answer"])      
    except:
        return "-"

def check_parentheses_and_keywords(row):
    if '(' in row and ')' in row and 'ukb' in row and 'field' in row:
        return row
    else:
        return "-"

 
 

def create_directory(directory):
    """Function to create a directory if it doesn't exist."""
    if not os.path.exists(directory):  # Checking if the directory doesn't exist
        os.makedirs(directory)  # Creating the directory if it doesn't exist
    return directory  # Returning the created or existing directory



def searchphenotypeandpopulation(phenotype,population):
    df = pd.read_csv("summary_statistics_table_export.tsv", sep='\t')
    savefile = phenotype
    phenotype = phenotype.replace("_"," ")

    df['reportedTrait'] = df['reportedTrait'].str.lower()
    df['efoTraits'] = df['efoTraits'].str.lower()
    df['discoverySampleAncestry'] = df['discoverySampleAncestry'].str.lower()
    df['UKB_FIELD'] = df['reportedTrait'].apply(check_parentheses_and_keywords)

    filter1 = df.copy()
    
    def fuzzy_similarity1(row):
        return fuzz.token_sort_ratio(phenotype, row['reportedTrait'])
    
    def fuzzy_similarity2(row):
        return fuzz.token_sort_ratio(population, row['discoverySampleAncestry'])
    
    filter1['FuzzySimilarity_'+phenotype] = filter1.apply(fuzzy_similarity1, axis=1)
    filter1 = filter1.sort_values(by='FuzzySimilarity_'+phenotype, ascending=False)
    phenotypethreshold = 50
    filter1 = filter1[filter1['FuzzySimilarity_'+phenotype]>phenotypethreshold]
    print(f"The number of GWAS after removing phenotypes with a fuzzy similarity score below {phenotypethreshold} is {len(filter1)}")
    
    if population is None:
        pass
    else:
        filter1['FuzzySimilarity_'+population] = filter1.apply(fuzzy_similarity2, axis=1)
        populationthreshold = 50
        filter1 = filter1[filter1['FuzzySimilarity_'+population]>populationthreshold]
        print(f"The number of GWAS after removing populations with a fuzzy similarity score below {populationthreshold} is {len(filter1)}")
        
    filter1["SearchPhenotype"] = phenotype



    # Turn of the following code if you do not care about the cases and controls.
    #question1 = "what is the number of controls mentioned in the text?"
    #filter1['CONTROLS'] = filter1["initialSampleDescription"].apply(getmeanswers_cases_controls,question = question1)

    #question1 = "what is the number of cases mentioned in the text?"
    #filter1['CASES'] = filter1["initialSampleDescription"].apply(getmeanswers_cases_controls,question = question1)

    #question1 = "what is the number of samples mentioned in the text?"
    #filter1['SAMPLES'] = filter1["initialSampleDescription"].apply(getmeanswers_samples,question = question1)


    filter1.to_csv(phenotype.replace(" ","_")+".csv")
    print(f"Total number of GWAS for Phenotype: {phenotype} and Population: {population} is {len(filter1)}")
    
    
    
    
    
    pass


import argparse

def main():
    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description='Search GWAS data based on phenotype and population.')

    # Add command-line arguments
    parser.add_argument('--phenotype', type=str, required=True, help='Specify the phenotype name.')
    parser.add_argument('--population', type=str,default=None, help='Specify the population.')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Access the values of the arguments
    phenotype = args.phenotype.lower()
    try:
        population = args.population.lower()
    except:
        population = args.population

    # Your processing logic based on the provided phenotype and population
    if phenotype:
        print(f"Phenotype: {phenotype}")
        print(f"Population: {population}")
        # Add your processing logic here
        
        searchphenotypeandpopulation(phenotype,population)

    else:
        print("Please provide phenotype")

if __name__ == "__main__":
    main()

exit(0)

 	 


