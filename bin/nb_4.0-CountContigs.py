
#################################################
### THIS FILE WAS AUTOGENERATED! DO NOT EDIT! ###
#################################################
# file to edit: dev_nb/4.0-CountContigs.ipynb
import json
import sys
import os
import os.path
from pathlib import Path
from subprocess import call
import subprocess
import sys # system libraries, like arguments (argv)
import re # regular expressions
import pandas as pd
import glob
import os.path
from pathlib import Path



from timeit import default_timer as timer
import time



def makeFolders(folder_list):
    for directory in folder_list:
        if not directory: # Make sure not an empty string ""
            continue
        else:
            if not os.path.exists(directory):
                os.makedirs(directory)

class fastQ():

    def __init__(self, NF_out, col_data):
        self.NF_out = NF_out
        self.wkdir = f'{self.NF_out}/unmapped/final'
        self.col_data = col_data
        self.df_col = pd.read_csv(self.col_data,sep="\t",names = ["sample","condition"])
        self.condition_list = self.df_col["condition"].unique()
        self.df_fq = pd.DataFrame()


    def getfqCountInitial(self, f):
        df = pd.read_csv(f,sep=" ",names=["count"])
        count = df["count"].iloc[0]
        count = int(int(count) / 4)
        name = os.path.basename(f)
        name = name.replace("_R1_fastqtrim_count.txt","")
        name = name.replace("_R1_fastq_count.txt","")
        return name, count



    def getfqCount(self, f):
        df = pd.read_csv(f,sep=" ",names=["count","name"])
        count = df["count"].iloc[0]
        count = int(int(count) / 4)
        name = os.path.basename(f)
        name = name.split("Unmapped")[0]
        return name, count



    def fastqtonumber(self):
        fq_path = f'{self.NF_out}/unmapped/final/filter/fastq/'
        start_fqs = sorted(glob.glob(f'{self.NF_out}/qc/fastq_count/*count.txt'))
        trim_fqs = sorted(glob.glob(f'{self.NF_out}/qc/fastqtrim_count/*count.txt'))
        unmapped_fqs = sorted(glob.glob(f'{fq_path}/initial/*count.txt'))
        afterbowtie_fqs = sorted(glob.glob(f'{fq_path}/bowtie2/*bowtie2_count.txt'))
        magicblast_fqs = sorted(glob.glob(f'{fq_path}/magicblast/*magicblast_count.txt'))

        df = pd.DataFrame(columns = ["sample","count","step"])
        f = unmapped_fqs[0]
        x = 0



        for f in start_fqs:
            name, count = self.getfqCountInitial(f)
            df.loc[x] = pd.Series( {"sample":name,"count":count, "step":"initial_reads"} )
            df["count"] = df["count"].astype(int)
            x = x + 1


        for f in trim_fqs:
            name, count = self.getfqCountInitial(f)
            df.loc[x] = pd.Series( {"sample":name,"count":count, "step":"trimmed_reads"} )
            df["count"] = df["count"].astype(int)
            x = x + 1



        for f in unmapped_fqs:
            name, count = self.getfqCount(f)
            df.loc[x] = pd.Series( {"sample":name,"count":count, "step":"unmapped_star"} )
            df["count"] = df["count"].astype(int)
            x = x + 1
        for f in afterbowtie_fqs:
            name, count = self.getfqCount(f)
            df.loc[x] = pd.Series( {"sample":name,"count":count, "step":"unmapped_bowtie2"} )
            x = x + 1
        for f in magicblast_fqs:
            name, count = self.getfqCount(f)
            df.loc[x] = pd.Series( {"sample":name,"count":count, "step":"unmapped_magicblast"} )
            x = x + 1

        df_star = df[df["step"]=="initial_reads"]
        df_merge = pd.merge(self.df_col,df_star, on="sample")
        df_merge = df_merge.rename(columns={"count":"initial_reads"})
        del df_merge["step"]

        df_temp = df[df["step"]=="trimmed_reads"]
        df_merge = pd.merge(df_merge,df_temp, on="sample")
        df_merge = df_merge.rename(columns={"count":"trimmed_reads"})
        del df_merge["step"]

        df_temp = df[df["step"]=="unmapped_star"]
        df_merge = pd.merge(df_merge,df_temp, on="sample")
        df_merge = df_merge.rename(columns={"count":"unmapped_star"})
        del df_merge["step"]

        df_temp = df[df["step"]=="unmapped_bowtie2"]
        df_merge = pd.merge(df_merge,df_temp, on="sample")
        df_merge = df_merge.rename(columns={"count":"unmapped_bowtie2"})
        del df_merge["step"]

        df_temp = df[df["step"]=="unmapped_magicblast"]
        df_merge = pd.merge(df_merge,df_temp, on="sample")
        df_merge = df_merge.rename(columns={"count":"unmapped_magicblast"})
        del df_merge["step"]
        df = df_merge


        df_initial = pd.read_csv(f"{self.NF_out}/qc/star_mapstats/starmapstats.txt",sep="\t")

        df = pd.merge(df,df_initial,on="sample")



        gb = df.groupby("condition").sum()
        gb["sample"] = gb.index
        gb["condition"] = gb.index

        df_all = pd.DataFrame(df.sum().iloc[2:])
        df_all = df_all.rename(columns={0:"all"}).T
        df_all["sample"] = "all"
        df_all["condition"] = "all"


        df = pd.concat([df_all,gb,df], sort=False)


        df["trimmed_reads_removed"] = df["initial_reads"] - df["trimmed_reads"]
        df["mapped_star"] = df["trimmed_reads"] - df["unmapped_star"]
        df["mapped_bowtie2"] = df["unmapped_star"] - df["unmapped_bowtie2"]

#         df["mapped_magicblast"] = df["unmapped_bowtie2"] - df["unmapped_magicblast"]
        df["unmapped_final"] = df["initial_reads"] - df["trimmed_reads_removed"] - df["mapped_star"] - df["mapped_bowtie2"]

        df["trimmed_reads_removed_perc"] = df["trimmed_reads_removed"] / df["initial_reads"]
        df["mapped_star_perc"] = df["mapped_star"] / df["initial_reads"]
        df["mapped_bowtie2_perc"] = df["mapped_bowtie2"] / df["initial_reads"]
#         df["mapped_magicblast_perc"] = df["mapped_magicblast"] / df["initial_reads"]
        df["unmapped_final_perc"] = df["unmapped_final"]  / df["initial_reads"]







        self.df_fq = df



        df = df[["sample","condition","initial_reads","trimmed_reads","trimmed_reads_removed",
                 "mapped_star",
                 "unmapped_star","mapped_bowtie2","unmapped_bowtie2","unmapped_final" ,
        "trimmed_reads_removed_perc",
        "mapped_star_perc",
        "mapped_bowtie2_perc" ,
        "unmapped_final_perc"   ]]

#         df["percent mapped star"] = df["mapped_star"] / df["input"]
#         df["percent mapped bowtie2"] = df["mapped_bowtie2"] / df["input"]
#         df["percent unmapped"] = df["unmapped_final"] / df["input"]


        df.to_csv(f'{fq_path}/fastq_amount_df.csv',sep="\t",index=None)



    def plotBarChartUnmapped(self):
        df = self.df_fq
        df = df[["sample","unmapped_star","unmapped_bowtie2","unmapped_magicblast"]]
        df.index = df["sample"]
        # del df_stacked["condition"]
        fig = df.plot(kind='bar',  figsize=(20, 16), fontsize=15).get_figure()

        self.barpath = f'{self.NF_out}/plots/barchart'
        makeFolders([f'{self.NF_out}/plots',self.barpath])
        fig.savefig(f'{self.barpath}/unmapped_reads.pdf')



#         ${unmappedPath}/final/filter/fastq/concatCondition"
#         wc -l group_${condition}_unmapped_R1.fastq > group_${condition}_unmapped_fastq_count.txt



class countContigs():
    def __init__(self, NF_out, col_data):
        self.NF_out = NF_out
        self.wkdir = f'{self.NF_out}/unmapped/final'
        self.col_data = col_data

    def countContig(self,f):
        df = pd.read_csv(f,sep="\t")
        length = len(df)
        return length

    def collectContigs(self):

#         out = glob.glob{f'{self.wkdir}/filter/contigs/
        levels = ["all","group","single"]
        removelist = ["9606","10090","Chordata","Viridiplantae","Artificial"]
        cols = ["sample","initial_contigs"]
        df = pd.DataFrame(columns=cols)


        samples = []
        lengths = []
        for level in levels:
            files = glob.glob(f"{self.wkdir}/filter/contigs/blast/initial_contigs/{level}/*_contigs.txt")

            for f in sorted(files):
                name = os.path.basename(f).replace("_unmapped_initial_contigs.txt","")
                length = self.countContig(f)
                samples.append(name)
                lengths.append(length)

        df["sample"] = samples
        df["initial_contigs"] = lengths

        for remove in removelist:
            lengths = []
            for level in levels:
                files = glob.glob(f"{self.wkdir}/filter/contigs/blast/Removed_afterBlast_{remove}/{level}/*.txt")
                for f in sorted(files):
                    name = os.path.basename(f).replace(f"_{remove}.txt","")

                    length = self.countContig(f)
                    lengths.append(length)
            if remove == "9606":
                remove = "human"
            if remove == "10090":
                remove = "mouse"
            df[remove] = lengths

        df.to_csv(f'{self.wkdir}/filter/contigs/blast/contigs_amount_df.txt',sep="\t",index=None)



class DarkBiomeCount():
    def __init__(self, NF_out):
        self.NF_out = NF_out
        self.wkdir = f'{self.NF_out}/unmapped/final'
        self.df = pd.DataFrame()
        self.df["file"] = ""
        self.df["initial"] = ""
        self.df["afterDust"] = ""
        self.df["afterVec"] = ""





    def readFiles(self, count_step):
        files = sorted(glob.glob(f"{self.wkdir}/darkbiome/{count_step}/*/*") )
        counts = []
        names = []
        for file in files:
            name = os.path.basename(file)
            x = self.fastaCount(file)
            counts.append(x)
            names.append(name)
        self.df[count_step] = counts
        self.df["file"] = names



    def fastaCount(self,file):
        with open(file,"r") as infile:
            x = 0
            for line in infile:
                if ">" in line:
                    x += 1
                    if "DUMMY" in line:
                        x -= 1
            return x

    def runDarkBiomeCount(self):
        self.readFiles(count_step="initial")
        self.readFiles(count_step="afterDust")
        self.readFiles(count_step="afterVec")
        self.df.to_csv(f"{self.wkdir}/darkbiome/contigs_amount_dark_df.txt",sep="\t", index=None)




class lastDBdiff():
    def __init__(self, NF_out):
        self.NF_out = NF_out
        self.wkdir = f'{self.NF_out}/unmapped/final/darkbiome'


    def dfDIFF(self,files,files2, outfile, level,trimend):
        for x in range(len(files)):
            f = files[x]
            f = os.path.basename(f)
            l = len(trimend)
            print("working on: ", f)
            name = f'{f[:-l]}'
            df = pd.read_csv(files[x],sep="\t")
            df2 = pd.read_csv(files2[x],sep="\t")
            L = len(df) - len(df2)
            outfile.write(f'{name}\t{level}\tBIOME\tRegular_vs_LASTDB\t{L}\n')



    def fastaCount(self,file):
        with open(file,"r") as infile:
            x = 0
            for line in infile:
                if ">" in line:
                    x += 1
                    if "DUMMY" in line:
                        x -= 1
            return x

    def fa_DIFF(self,files,files2, outfile, level,trimend):
        for x in range(len(files)):

            f = files[x]
            f = os.path.basename(f)
            l = len(trimend)
            name = f'{f[:-l]}'
            fa_count = self.fastaCount(files[x])


#             df2 = pd.read_csv(files2[x],sep="\t")
#             count2 = len(df2)
            fa_count2 = self.fastaCount(files2[x])
            L = fa_count - fa_count2
            outfile.write(f'{name}\t{level}\t{fa_count}\t{fa_count2}\t{L}\n')

    def amountDarkLASTDB(self):
        print("Running ")
          #Same for Dark genome

        f_single = sorted(glob.glob(f'{self.wkdir}/afterVec/single/*biome.fasta'))
        single_FG = sorted(glob.glob(f'{self.wkdir}/lastdb/lastdb_fasta/group_single/single/*.fasta'))

        f_group = sorted(glob.glob(f'{self.wkdir}/afterVec/group/*biome.fasta'))
        group_FA = sorted(glob.glob(f'{self.wkdir}/lastdb/lastdb_fasta/all_group/group/*.fasta'))

        f_all = sorted(glob.glob(f'{self.wkdir}/afterVec/all/*biome.fasta'))
        all_FG = sorted(glob.glob(f'{self.wkdir}/lastdb/lastdb_fasta/all_group/all/*.fasta'))


        fout = f'{self.wkdir}/contigs_amount_dark_lastdb_df.txt'
        with open(fout,"w") as outfile:

            trimend = "_darkbiome.fasta"
            outfile.write(f'file\tlevel\t\tamount_afterDust\tamount_lastdb\tamount_removed\n')
            self.fa_DIFF(f_single,single_FG,outfile,"single_FG",trimend)
            self.fa_DIFF(f_group,group_FA,outfile,"group_FA",trimend)
            self.fa_DIFF(f_all,all_FG,outfile,"all_FG",trimend)



class amountDark():
    def __init__(self, NF_out):
        self.NF_out = NF_out
        self.wkdir = f'{self.NF_out}/unmapped/final/darkbiome'


    def readDF(self, f):
        try:
            df = pd.read_csv(f,sep="\t")
            df = df[~df["contig"].str.contains("DUMMY")]
#             DUMMYDARK

            count_all = len(df)
        except:
            count_all = 0
        return count_all



    def runDark(self):

        x_nohit = self.readDF(f"{self.wkdir}/lastdb/FINAL_Nohits/all_group/all/all_darkbiome_df.txt" )

        tophits = "lastdb/FINAL_Tophits/all_group/all/blastXhits"
        x_art = self.readDF(f"{self.wkdir}/{tophits}/all_darkb_Artificial.txt" )
        x_cho = self.readDF(f"{self.wkdir}/{tophits}/all_darkb_Chordata.txt" )
        x_vir = self.readDF(f"{self.wkdir}/{tophits}/all_darkb_Viridiplantae.txt" )
        x_top = self.readDF(f"{self.wkdir}/lastdb/FINAL_Tophits/all_group/all/all_darkbiome_df.txt")

        with open(f"{self.wkdir}/contigs_amount_dark_afterblast.txt","w") as outfile:
            outfile.write(f"step\tamount\nFinal_No_Hits\t{x_nohit}\nFinal_Artificial\t{x_art}\nFinal_Chordata\t{x_cho}\nFinal_Viridiplantae\t{x_vir}\nFinal_Hits\t{x_top}"
            )



import argparse
def parse_arguments():
        parser = argparse.ArgumentParser(description='filter quantify and graph')
        parser.add_argument('--col_data', action= 'store', metavar='col_data')
        parser.add_argument('--Nextflow_Out', action= 'store', metavar='--Nextflow_Out')
        args = parser.parse_args()
        return args

if __name__=="__main__":
    args = parse_arguments()
    col_data = args.col_data
    NF_out = args.Nextflow_Out


    F = fastQ(NF_out = NF_out, col_data = col_data)
    F.fastqtonumber()


    C = countContigs(NF_out = NF_out, col_data = col_data)
    C.collectContigs()


    D = DarkBiomeCount(NF_out=NF_out)
    D.runDarkBiomeCount()

    D = lastDBdiff(NF_out=NF_out)
    D.amountDarkLASTDB()

    D = amountDark(NF_out=NF_out)
    D.runDark()