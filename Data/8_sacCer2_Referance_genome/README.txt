
This directory contains the June 2008 assembly of the S. cerevisiae genome
(sacCer2, SGD June 2008) in one gzip-compressed FASTA file per chromosome.

The data is based on sequence dated June 2008 in the Saccharomyces Genome
Database (http://www.yeastgenome.org/) and was obtained from the site
http://downloads.yeastgenome.org/sequence/genomic_sequence/chromosomes/fasta/
The S288C strain was used in this sequencing project.

Files included in this directory:

  - chr*.fa.gz, 2micron.fa.gz: compressed FASTA sequence of each chromosome.
	No masking has been applied to this sequence.
  - md5sum.txt: md5sum checksums for the files in this directory.
------------------------------------------------------------------
If you plan to download a large file or multiple files from this 
directory, we recommend that you use ftp rather than downloading the 
files via our website. To do so, ftp to hgdownload.cse.ucsc.edu, then 
go to the directory goldenPath/sacCer2/chromosomes. To download multiple 
files, use the "mget" command:

    mget <filename1> <filename2> ...
    - or -
    mget -a (to download all the files in the directory)

Alternate methods to ftp access.
    
Using an rsync command to download the entire directory:
    rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/sacCer2/chromosomes/ .
For a single file, e.g. chrM.fa.gz
    rsync -avzP 
        rsync://hgdownload.cse.ucsc.edu/goldenPath/sacCer2/chromosomes/chrM.fa.gz .
    
Or with wget, all files:
    wget --timestamping 
        'ftp://hgdownload.cse.ucsc.edu/goldenPath/sacCer2/chromosomes/*'
With wget, a single file:
    wget --timestamping 
        'ftp://hgdownload.cse.ucsc.edu/goldenPath/sacCer2/chromosomes/chrM.fa.gz' 
        -O chrM.fa.gz
    
To uncompress the fa.gz files:
    gunzip <file>.fa.gz


All the tables in this directory are freely available for public use.
