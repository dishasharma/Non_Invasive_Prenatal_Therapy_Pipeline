import os
import pandas as pd
import numpy as np
import glob
#mainbeddir = 'bedfiles/'
#bed_dir_list = os.listdir(mainbeddir)
#bed_listing = glob.glob(os.path.join(mainbeddir, '*.bed'))
#bed_listing.sort()
gc_con_arr = []
fname = "binsize_method/genome_20kb.bed"
#for fname in bed_listing:
cmnd = "bedtools nuc -fi hg38.fa -bed " + fname + " | awk \'{ print $1\":\"$2\"-\"$3\"\t\"$5}\'"
gc_content = os.popen(cmnd)
#gc_con_arr.append(fname + gc_content)
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
print(dff_l[0:5])