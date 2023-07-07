import os
import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt

#GC correct the bin content


#open both bed file and fasta here
mainbeddir = 'bam_files/'                                                  #functional directory here
#mainfastadir = 'fastafiles/'
bed_dir_list = os.listdir(mainbeddir)
#fasta_dir_list = os.listdir(mainfastadir)
#xdir = 'xyzabc.kmerlists'
#os.makedirs(xdir)
bed_listing = glob.glob(os.path.join(mainbeddir, '*bincount_50kb.bed'))               #creates listing of bed files in directory
#fasta_listing = glob.glob(os.path.join(mainfastadir, '*.fasta'))
bed_listing.sort()
#fasta_listing.sort()
read_list = []
bin_list = []
gc_list = []

cfname = "binsize_method/genome_50kb.bed"

cmnd = "bedtools nuc -fi hg38.fa -bed " + cfname + " | awk \'{ print $1\":\"$2\"-\"$3\"\t\"$5}\'"
gc_content = os.popen(cmnd)

x = gc_content.read()
b_list = x.strip().split()
length = len(b_list)
print(length)
x_list = []
loclen = int(length/2)
for i in range(0,loclen):
    abc = []
    abc.append(b_list[2*i])
    abc.append(b_list[(2*i)+1])
    x_list.append(abc)

dff_f = pd.DataFrame(data = x_list[1:],columns = x_list[0])
dff_l = dff_f.values.tolist()
for b in dff_l:
    gc_list.append(b[1])

for fname in bed_listing:
    prefix = os.path.splitext(os.path.split(fname)[-1])[0]          #prefix maker
    #list_fname = os.path.join(xdir, prefix + '_list.txt')
    content = []
    content_chrname = []
    content_readcount = []
    with open(fname)as f:
        for line in f:
            content.append(line.strip().split())
    for x in content:
        chrname = x[0]+":"+x[1]+"-"+x[2]
        content_chrname.append(chrname)
        content_readcount.append(int(x[3]))
    read_list.append(content_readcount)

bin_list = content_chrname

x_df = pd.DataFrame(data = read_list, index = list(range(1,(len(bed_listing)+1))), columns = bin_list)

print(pd.Series(data = gc_list))
#determine number of reads per bin

plt.scatter(gc_list,read_list[0])
plt.ylim(top = 500, bottom = 0)
plt.show()
#get GC content of each bin

xer = 1
#LOESS plot on the no. of reads vs GC content
for i in read_list:
    lowess = sm.nonparametric.lowess(i, gc_list, frac=.3)
    lowess_x = list(zip(*lowess))[0]
    lowess_y = list(zip(*lowess))[1]

# run scipy's interpolation. There is also extrapolation I believe
    f = interp1d(lowess_x, lowess_y, bounds_error=False)

    xnew = [b/10. for b in range(400)]

# this this generate y values for our xvalues by our interpolator
# it will MISS values outsite of the x window (less than 3, greater than 33)
# There might be a better approach, but you can run a for loop
#and if the value is out of the range, use f(min(lowess_x)) or f(max(lowess_x))
    ynew = f(xnew)

    plt.figure(xer)
    plt.plot(gc_list, i, 'o')
    plt.plot(lowess_x, lowess_y, '*')
    plt.plot(xnew, ynew, '-')
    xer = xer + 1

plt.show()
#determine median read count of all bins (tbd if avg or median)


#GC correct the read counts


#output to file