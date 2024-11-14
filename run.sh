#Download and execute Fastqc
wget 'https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip'
unzip -d . fastqc_v0.11.8.zip
cd Fastqc
chmod +x fastqc
find ~/fastqfolder -name *.fastq* | parallel "fastqc {}"

#Quality trimming of reads
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip -d . Trimmomatic-0.39.zip
cd Trimmomatic-0.39
find ~/fastqfolder -name *.fastq* | parallel "*** {}" #Input single-end/paired-end reads then run command accordingly


#Download the Reference Genome from UCSC genome browser
wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'   
gunzip hg38.fa.gz

#Install BWA and Index reference genome
git clone http://github.com/lh3/bwa.git
cd bwa
make

#Index the Reference genome
bwa index -p hg38ref_index -a bwtsw hg38.fa

#Alignment
./bwa mem hg38ref_index *.fastq.gz *.sam #Read Fastq file here and output file accordingly

#Aignment post-processing
git clone http://github.com/samtools/samtools.git
./configure
make
make install

#Sam to bam
samtools view -S *.sam -b -o *.bam
samtools sort *.bam -o *.sorted
samtools index *.sorted

###########Count-based method to identify chromosomal aneuploidy

#Python code to extract per chromosome reads
####################### ADD PYTHON ENVIRONMENT HERE ###############
# importing the pysam package that enables the parsing of BAM files to python
import pysam
#built-in method for calling the alignment file
samfile=pysam.AlignmentFile("*")
# searches for the specified strings in the alignment files
list=(“chr1”,”chr2”,”chr3”,”chr4”,”chr5”,”chr6”,”chr7”,”chr8”,”chr9”,”chr10”,”chr11
”,”chr12”,”chr13”,”chr14”,”chr15”,”chr16”,"chr17","chr18","chr19","chr20",”chr21”,
”chr22”,”chrX”)
# employs the fetch method in a loop to call forth the unique and specific read
mapping to a particular chromosome
for i in list:
iter=samfile.fetch(i)
print (str(i))
count=0
for x in iter:
print (str(x))
count+=1
print (count)
#end

###########Size-based method to identify chromosomal abnormality
#Get chromosome sizes for hg38 reference
wget –timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes'

#Binning the reference genome
git clone https://github.com/arq5x/bedtools2.git
./configure
make
bedtools makewindows -g hg38.chrom.sizes -w *size > genome_*.bed

#BamtoBed
bamToBed -i *.sorted > *.bed

#Genome coverage per sample
bedtools coverage -a genome_*.bed -b *.bed –counts > *.bed

#Awk installation
git clone https://github.com/onetrueawk/awk.git
make
awk ‘{print $1"\t"$2"\t"$3"\t"$4"\t"$1":"$2"-"$3"\t"FILENAME }’ *.bed > *.bed

cat *.bed > merged_*.bed


#Normalize by Chromosome size and Read Count (divide by chromosomal size and read count)


###########K-mer method to identify chromosomal abnormality
mkdir Reference_genome_split
cd Reference_genome_split
awk '{ if(NR==1 && $0 ~/>/){c=$0;gsub(">","",c);print $0 > c".fa";} else if($0~/>/ && NR>1){close(c".fa"); c=$0;gsub(">","",c);print $0 > c".fa"} else { print $0 > c".fa";}}' ../hg38.fa

#install Biopython 

pip install biopython


#Make K-mers for reference genome
def fasta_reader(filename):
    from Bio.SeqIO.FastaIO import FastaIterator
    with open(filename) as handle:
        for record in FastaIterator(handle):
            yield record


maindir = 'xyzabc'                                                  #functional directory here
dir_list = os.listdir(maindir)
#xdir = 'xyzabc.kmerlists'
#os.makedirs(xdir)
listing = glob.glob(os.path.join(maindir, '*.fasta'))               #creates listing of fasta files in directory


for fname in listing:
    prefix = os.path.splitext(os.path.split(fname)[-1])[0]          #prefix maker
    list_fname = os.path.join(xdir, prefix + '_referenceKmer.txt')           #text file namer
    for entry in fasta_reader(fname):
        print(str(entry.id))                                        #This is header of fasta entry
        print(str(entry.seq))                                       #This is sequence of specific fasta entry
        seq_list = str(entry.seq)
    kmer_list = list(dict(make_kmer_list(10, seq_list)))            #k-mer maker with list-dict conversion to remove redundancy
    f_h = open(list_fname,'w+')
    for item in kmer_list:
        f_h.write("{}\n".format(item))                              #writing to text file
    f_h.close()
    


def make_kmer_list (k, alphabet):

    # Base case.
    if (k == 1):
        return(alphabet)

    # Handle k=0 from user.
    if (k == 0):
        return([])

    # Error case.
    if (k < 1):
        sys.stderr.write("Invalid k=%d" % k)
        sys.exit(1)

    # Precompute alphabet length for speed.
    alphabet_length = len(alphabet)

    # Recursive call.
    return_value = []
    for kmer in make_kmer_list(k-1, alphabet):
        for i_letter in range(0, alphabet_length):
            return_value.append(kmer + alphabet[i_letter])
              
    return(return_value)

#Substract the non-unique ones


#Make K-mer for fastq files
import Bio.SeqIO
import os
#Make K-mers for reference genome
def fastq_reader(filename):
    for seq_record in Bio.SeqIO.parse(filename, "fastq"):
        seq = seq_record
    yield seq


maindir = 'xyzabc'                                                  #functional directory here
dir_list = os.listdir(maindir)
#xdir = 'xyzabc.kmerlists'
#os.makedirs(xdir)
listing = glob.glob(os.path.join(maindir, '*.fastq'))               #creates listing of fastq files in directory


for fname in listing:
    prefix = os.path.splitext(os.path.split(fname)[-1])[0]          #prefix maker
    list_fname = os.path.join(xdir, prefix + '_sampleKmer.txt')           #text file namer
    for entry in fastq_reader(fname):
        print(str(entry.id))                                        #This is header of fastq entry
        print(str(entry.seq))                                       #This is sequence of specific fastq entry
        seq_list = str(entry.seq)
    kmer_list = list(dict(make_kmer_list(10, seq_list)))            #k-mer maker with list-dict conversion to remove redundancy
    f_h = open(list_fname,'w+')
    for item in kmer_list:
        f_h.write("{}\n".format(item))                              #writing to text file
    f_h.close()
    


def make_kmer_list (k, alphabet):

    # Base case.
    if (k == 1):
        return(alphabet)

    # Handle k=0 from user.
    if (k == 0):
        return([])

    # Error case.
    if (k < 1):
        sys.stderr.write("Invalid k=%d" % k)
        sys.exit(1)

    # Precompute alphabet length for speed.
    alphabet_length = len(alphabet)

    # Recursive call.
    return_value = []
    for kmer in make_kmer_list(k-1, alphabet):
        for i_letter in range(0, alphabet_length):
            return_value.append(kmer + alphabet[i_letter])
              
    return(return_value)


#Counting per chromosome K-mers for each sample
import os
import glob

row_count = 0
maindir = ""
output = []
#sample file handler
sample_file_list = glob.glob(os.path.join(maindir, '*_sampleKmer.txt'))
sample_file_list.sort()
#reference file handler
ref_file_list = glob.glob(os.path.join(maindir, '*_referenceKmer.txt'))
ref_file_list.sort()
for file in sample_file_list:
    col_count = 0
    out_1_row = []
    for ref in ref_file_list:
        print(file)
        print(ref)
        cmd = "grep -c -f " + file + " " + ref         #command line feeder
        p = os.popen(cmd)                              #output as string
        #print(type(p.read()))                         #type verification
        count_temp = int(p.read())                     #integer conversion
        #print(type(count_temp))                       #re type check
        out_1_row.append(count_temp)                   #column append
        col_count = col_count + 1                      #column count increment
    output.append(out_1_row)                           #row append
    row_count = row_count + 1                          #row count increment
	
print(output)


