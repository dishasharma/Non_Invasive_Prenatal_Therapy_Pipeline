import Bio.SeqIO
import os
import glob
#Make K-mers for reference genome
def fastq_reader(filename):
    for seq_record in Bio.SeqIO.parse(filename, "fasta"):
        seq = seq_record
    yield seq

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
    kmer_count = alphabet_length - k + 1

    # Recursive call.
    return_value = []
    for i in list(range(0,kmer_count)):
    	temp = i+k
    	return_value.append(alphabet[i:temp])
              
    return(return_value)

maindir = 'xyzabc'                                                  #functional directory here
dir_list = os.listdir(maindir)
#xdir = 'xyzabc.kmerlists'
#os.makedirs(xdir)
listing = glob.glob(os.path.join(maindir, '*.fasta'))               #creates listing of fastq files in directory


for fname in listing:
    prefix = os.path.splitext(os.path.split(fname)[-1])[0]          #prefix maker
    list_fname = os.path.join(maindir, prefix + '_list.txt')           #text file namer
    for entry in fastq_reader(fname):
        print(str(entry.id))                                        #This is header of fastq entry
        print(str(entry.seq))                                       #This is sequence of specific fastq entry
        seq_list = str(entry.seq)
    kmer_list = list(dict.fromkeys(make_kmer_list(10, seq_list)))            #k-mer maker with list-dict conversion to remove redundancy
    f_h = open(list_fname,'w+')
    for item in kmer_list:
        f_h.write("{}\n".format(item))                              #writing to text file
    f_h.close()
    


