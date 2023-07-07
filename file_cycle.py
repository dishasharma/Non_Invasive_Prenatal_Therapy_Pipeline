import os
import glob
#file cycler
maindir = 'bam_files/'                                                  #functional directory here
dir_list = os.listdir(maindir)
xdir = "split_bam_files"
os.system("mkdir "+xdir)
listing = glob.glob(os.path.join(maindir, '*bincount_50kb.bed'))               #creates listing of bed files in directory
listing.sort()


for fname in listing:
    prefix = os.path.splitext(os.path.split(fname)[-1])[0]          #prefix maker
    list_fname = os.path.join(xdir, prefix[:9])
    cmnd = "awk \'{print $0 >> \"" + list_fname +"_\"$1\".bed\"}\' " + fname
    print(cmnd)
    os.system(cmnd)