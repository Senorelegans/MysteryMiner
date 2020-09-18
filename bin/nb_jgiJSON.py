
#################################################
### THIS FILE WAS AUTOGENERATED! DO NOT EDIT! ###
#################################################
# file to edit: dev_nb/jgiJSON.ipynb
import json
import sys
import os
import os.path
from pathlib import Path
from subprocess import call
import sys # system libraries, like arguments (argv)
import re # regular expressions
import pandas as pd


from timeit import default_timer as timer
import time


def removeHuman(f1,fout):
    my_cols=[
    "contig",
    "blast_pident",
    "blast_sseqid",
    "blast_evalue",
    "taxid",
    "blast_sscinames",
    "blast_scomnames",
    "blast_sskingdoms",
    "blast_stitle"
    ]
    df=pd.read_csv(f1, dtype={
                            "contig":str,
                            "blast_pident":str,
                            "blast_sseqid":str,
                            "blast_evalue":str,
                            "taxid":str,
                            "blast_sscinames":str,
                            "blast_scomnames":str,
                            "blast_sskingdoms":str,
                            "blast_stitle":str
                            }, sep='\t', names=my_cols)

    del df["blast_sscinames"]
    del df["blast_scomnames"]
    del df["blast_sskingdoms"]
    df["filename"] = fout

    df['contig'] = df.drop_duplicates(subset=['contig'], keep="first")
    df.to_csv(f'{fout}_initial_contigs.txt',index=None,sep="\t")


    def getDiff(taxid,df):
        df2 = df[df["taxid"]== taxid]
        df2.to_csv(f'{fout}_{taxid}.txt',index=None,sep="\t")
        df = df[df["taxid"]!= taxid]
        return df

    df = getDiff(taxid="9606",df=df) # Remove human
    df = getDiff(taxid="10090",df=df) # Remove mouse
    return df

def firstJGIQuery(df,fout):

    taxids=df['taxid'].unique()
    df_ids = pd.DataFrame(taxids)
    df_ids = df_ids.rename(columns = {0:"taxid"})

    df_ids["jgi_tab"] = "NA"
    df_ids["jgi_json"] = "NA"
    i = 0
    R = len(df_ids)
    for x in range(R):
        print(f'Querying {i+1} of {R} taxids')
        ids = df_ids["taxid"].iloc[x]
        Q1 = "curl https://taxonomy.jgi-psf.org/sc/id/" + str(ids)
        result = os.popen(Q1).read()
        df_ids["jgi_tab"].iloc[x] = result

        #Get json return

        jasonQ1 = "curl https://taxonomy.jgi-psf.org/id/" + str(ids)
        result = os.popen(jasonQ1).read()
        j = json.loads(result)
        j = j[list(j.keys())[0]]
        j = str(j)
        df_ids["jgi_json"].iloc[x] = j

        i+= 1



    df = pd.merge(df,df_ids,on="taxid")
    # To remove "Chordata","Viridiplantae",
    # Will remove at jgi json step "synthetic","artificial","PREDICTED"
    df_chordata = df[df["jgi_tab"].str.contains("Chordata")==True]
    df_virid = df[df["jgi_tab"].str.contains("Viridiplantae")==True]

    df_chordata.to_csv(f'{fout}_Chordata.txt',sep="\t",index=None)
    df_virid.to_csv(f'{fout}_Viridiplantae.txt',sep="\t",index=None)

    df_artificial = df[df["jgi_tab"].str.contains("synth|vector|Vector|artificial")==True]
    df_artificial.to_csv(f'{fout}_Artificial.txt',sep="\t",index=None)

    df = df[df["jgi_tab"].str.contains("Chordata|Viridiplantae")==False]
    df = df[df["jgi_tab"].str.contains("synth|vector|Vector|artificial")==False]

    return df


def matchFasta(df,fasta,fout):
    c11 = fout+"_newfasta1.fast"
    c12 = fout+"_newfasta2.fast"
    with open(fasta,'r', newline="\n") as infile, open(c11, 'w') as outfile:
        for line in infile:
            if ">" in line:
                line = line.replace('\n', '\t')
                line = line.replace('>', '\n')
                outfile.write(line)
            else:
                line = line.replace('\n', 'RTN666')
                outfile.write(line)
        infile.close()
        outfile.close()

    with open(c11,'r', newline="\n") as infile, open(c12, 'w') as outfile:
        for line in infile:
            if "NODE" not in line:
                line = line.replace('\n', '')
                line = line.replace('\r', '')
                outfile.write(line)
            else:
                outfile.write(line)
        infile.close()
        outfile.close()

    my_cols=['contig', 'fasta_with_returnsAdded']
    df_fa=pd.read_csv(c12, index_col=['contig'], sep='\t', names=my_cols)
    df = pd.merge(df, df_fa, on='contig') # Dont run this cell more than once
    df["fasta"] = df['fasta_with_returnsAdded'].str.replace("RTN666","")
    df["length"] = df["fasta"].str.len()
    df.to_csv(fout+"_jgi_filtered_df.txt",sep="\t",index=None)

    # Write out fasta
    with open(fout+"_jgi_filtered.fasta","w") as outfile:
        for x in range(len(df)):
            nodename = ">"+df["contig"].iloc[x]+"\n"
            fa = df['fasta_with_returnsAdded'].iloc[x]
            fa = fa.replace('RTN666', '\n')
            outfile.write(nodename+fa)


    return df


import ast

def toJSON(df,fout):
    df_J = df[["taxid","jgi_json"]]
    del df["jgi_json"]

    df.to_json(f"{fout}_temp.json",orient='index')
    with open(f"{fout}_temp.json") as fo:
        contig_d = json.load(fo)
        for c in contig_d.keys():
            d = contig_d[c]
            taxid = d["taxid"]
            df_temp = df_J[df_J["taxid"]==taxid]
            j = df_temp["jgi_json"].iloc[0]
            jdic = ast.literal_eval(j)
            d["jgi_json"] = jdic

        with open(f"{fout}_jgi.json", "w") as write_file:
            json.dump(contig_d, write_file)


import argparse
def parse_arguments():
        parser = argparse.ArgumentParser(description='filter fastas and turn to json')
        parser.add_argument('--file', action= 'store', metavar='file')
        parser.add_argument('--fasta', action= 'store', metavar='fasta')
        args = parser.parse_args()
        return args


if __name__=="__main__":
    args = parse_arguments()
    f = args.file
    fasta = args.fasta

    #Run all the commands
    fbase = os.path.basename(f)
    fout = fbase[:-4]
    df = removeHuman(f,fout)
    df = firstJGIQuery(df,fout)
    df = matchFasta(df,fasta,fout)
    toJSON(df,fout)

