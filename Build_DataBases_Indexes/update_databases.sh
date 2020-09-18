#!/usr/bin/env bash



# Make folders to put databases in
mkdir -p db
mkdir -p db/human_genomic
mkdir -p db/nt
mkdir -p db/nr


# Go to database top folder
cd db



#human nt database
cd human_genomic
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/human_genomic**tar.gz
for f in *.tar.gz; do tar -xvzf "$f"; done
cd ..


#all-organisms nt database
cd nt
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt**tar.gz
for f in *.tar.gz; do tar -xvzf "$f"; done
cd ..


#all-organisms nr database for blastx
cd nr
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr**tar.gz
for f in *.tar.gz; do tar -xvzf "$f"; done


echo done downloading databases
