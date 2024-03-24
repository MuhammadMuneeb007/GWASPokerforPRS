import pandas as pd
import numpy as np
import sys
import os
import argparse


def downloader(downloadurl,savingdirec):
    os.system("wget  "+"\"" +downloadurl+"\""+ " -P "+savingdirec)

    

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Download the full GWAS files!')

    parser.add_argument('--processedfile', type=str, help='Download the full GWAS files!')
    parser.add_argument('--indexer', type=int, help='Download the full GWAS files!')
    

    args = parser.parse_args()
    gwas_file = args.processedfile
    indexer = args.indexer

    data = pd.read_csv(gwas_file,encoding='utf-8') 
    print(data.head())
    
    names = data["Name"].values
    downloadlink = data["Download Link"].values
    print(names)
    print(downloadlink)
    print("Downloading ",len(pd.read_csv(gwas_file))," GWAS files from GWAS catalog!")
    
    #for loop in range(0,len(names)):
    import sys
    
    indexer = int(indexer)-1

    print("Downloading GWAS file for ",names[indexer]," from ",downloadlink[indexer])
    savingdirec = names[indexer] 
    downloadurl = downloadlink[indexer]
    downloader(downloadurl,savingdirec)
    #exit(0)










        



