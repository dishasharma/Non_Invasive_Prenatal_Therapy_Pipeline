import os
import glob
import pandas as pd

chr_dir = "binsize_method/test"
chr_listing = glob.glob(os.path.join(chr_dir, 'chr*.bed'))
chr_listing.sort()
print(pd.Series(data = chr_listing))

pre_list = []
ful_dir = "bam_files"
ful_listing = glob.glob(os.path.join(ful_dir, '*bincount_50kb.bed'))
for i in ful_listing:
    pre_ful = os.path.splitext(os.path.split(i)[-1])[0]          #prefix maker
    pre_list.append(pre_ful[:9])
pre_list.sort()

sam_dir = "split_bam_files"
sample_list = []
for s in pre_list:
    regex_file = s + '*.bed'
    sample_x = glob.glob(os.path.join(sam_dir, regex_file))
    sample_x.sort()
    sample_list.append(sample_x)

print(pd.Series(data = sample_list[0]))