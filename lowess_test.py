from math import ceil
import numpy as np
from scipy import linalg
import itertools
import pandas as pd
import os
import glob
import collections
import pylab as pl
import re


def lowess(x, y, f=2. / 3., iter=3):
    """lowess(x, y, f=2./3., iter=3) -> yest
    Lowess smoother: Robust locally weighted regression.
    The lowess function fits a nonparametric regression curve to a scatterplot.
    The arrays x and y contain an equal number of elements; each pair
    (x[i], y[i]) defines a data point in the scatterplot. The function returns
    the estimated (smooth) values of y.
    The smoothing span is given by f. A larger value for f will result in a
    smoother curve. The number of robustifying iterations is given by iter. The
    function will run faster with a smaller number of iterations.
    """
    n = len(x)
    r = int(ceil(f * n))
    h = [np.sort(np.abs(x - x[i]))[r] for i in range(n)]
    w = np.clip(np.abs((x[:, None] - x[None, :]) / h), 0.0, 1.0)
    w = (1 - w ** 3) ** 3
    yest = np.zeros(n)
    delta = np.ones(n)
    for iteration in range(iter):
        for i in range(n):
            weights = delta * w[:, i]
            b = np.array([np.sum(weights * y), np.sum(weights * y * x)])
            A = np.array([[np.sum(weights), np.sum(weights * x)],
                          [np.sum(weights * x), np.sum(weights * x * x)]])
            beta = linalg.solve(A, b)
            yest[i] = beta[0] + beta[1] * x[i]

        residuals = y - yest
        s = np.median(np.abs(residuals))
        delta = np.clip(residuals / (6.0 * s), -1, 1)
        delta = (1 - delta ** 2) ** 2

    return yest

#BEGINNING OF GC EXTRACTION
chr_dir = "binsize_method/test"
chr_listing = glob.glob(os.path.join(chr_dir, 'chr*.bed'))
chr_listing.sort()
final_gc = []
for cfname in chr_listing:
    gc_list = []
    cmnd = "bedtools nuc -fi hg38.fa -bed " + cfname + " | awk \'{ print $1\":\"$2\"-\"$3\"\t\"$5}\'"
    gc_content = os.popen(cmnd)

    x = gc_content.read()
    b_list = x.strip().split()
    length = len(b_list)
    #print(length)
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
        gc_list.append(float(b[1]))
    final_gc.append(gc_list)
#END OF GC EXTRACTION

#BEGINNING OF PREFIX TABULATION
pre_list = []
ful_dir = "bam_files"
ful_listing = glob.glob(os.path.join(ful_dir, '*bincount_50kb.bed'))
for i in ful_listing:
    pre_ful = os.path.splitext(os.path.split(i)[-1])[0]          #prefix maker
    pre_list.append(pre_ful[:9])
pre_list.sort()
#END OF PREFIX TABULATION

rex_sample = "^.*con[0-9]{1,2}.*chr[0-9XY]{1,2}.bed$"
rex_v = re.compile(rex_sample)
#BEGINNING OF BED FILE PROCESSING
sam_dir = "split_bam_files"
sample_list = []
for s in pre_list:
    regex_file = s + '*.bed'
    sample_x = glob.glob(os.path.join(sam_dir, regex_file))
    sample_x.sort()
    sample_list.append(sample_x)

#print(sample_list)
#fname = "binsize_method/test/chr1con1-bincount.bed"
#split into either chr_section wise or sample wise
sample_list_2 = [e for e in sample_list if e]

read_list = []

corr_final = []

corr_bin_fin = []
for s_x in sample_list_2:
    i_gc = 0
    corr_r_l = []
    corr_bin_l = []
    #print(s_x[0])
    for s_y in s_x:
        if rex_v.match(s_y):
            print(s_y)

            content = []
            content_binname = []
            content_readcount = []
            with open(s_y)as f:
                for line in f:
                    content.append(line.strip().split())
            for x in content:
                binname = x[0]+":"+x[1]+"-"+x[2]
                content_binname.append(binname)
                content_readcount.append(int(x[3]))
            #read_list.append(content_binname)

#END OF BED FILE PROCESSING

#BEGINNING OF DATA ARRANGEMENT AND LOESS PLOTTING
            import math
            x_lu = final_gc[i_gc]
            y_lu = content_readcount
            print(len(x_lu))
            print(len(y_lu))
            
            x_l = []
            y_l = []
            if len(x_lu) == len(y_lu):
                print("true")
            x_wit = []
            for i in range(0,len(x_lu)):
                if x_lu[i] != 0 and y_lu[i] != 0:
                    x_l.append(x_lu[i])
                    y_l.append(y_lu[i])
                    x_wit.append(content_binname[i])
            #print(x_l)
            #print(y_l)
            #print(len(x_lu))
            
            n = len(x_l)
            x = np.asarray(x_l)
            y = np.asarray(y_l)
            #print(type(y))
            ymean = np.mean(y_l)
            y_m = [ymean]*n
            f = 0.25
            yest = lowess(x, y, f=f, iter=3)
            corr_f = []
            for i in yest:
                c_f = i - ymean
                corr_f.append(c_f)

            y_p = []
            for j in range(0,n):
                c_y = y[j] - corr_f[j]
                y_p.append(c_y)
            x_s = []
            yest_s = []
            dic_val = dict(zip(x,yest))

            for key, val in sorted(dic_val.items()):
                x_s.append(key)
                yest_s.append(val)

            #pl.clf()
            #pl.scatter(x, y, label='y noisy', color = 'red')
            #pl.scatter(x, y_p, label='y correct', color = 'black')
            #pl.plot(x_s, yest_s, label='y pred')
            #pl.plot(x, y_m, label='y mean', color = 'green')
            #pl.ylim(0,200)
            #pl.legend()
            #pl.show()
#END OF DATA ARRANGEMENT AND LOESS PLOTTING

#BEGINNING OF DATA FORMATTING FOR PASS ON TO chi-sq-VR
            '''
            read_corr = []
            ye_c = 0
            for i in range(len(x_lu)):
                if x_lu[i] != 0 and y_lu[i] != 0:
                    read_corr.append(y_p[ye_c])
                    ye_c = ye_c + 1
                else:
                    read_corr.append(x_lu[i])
            #print(len(read_corr))
            '''
            corr_r_l.append(y_p)
            corr_bin_l.append(x_wit)
        i_gc = i_gc + 1
    flattened_list  = [y for x in corr_r_l for y in x]
    flattened_list_bin = [y for x in corr_bin_l for y in x]
    corr_final.append(flattened_list)
    corr_bin_fin.append(flattened_list_bin)

corr_final_2 = [e for e in corr_final if e]
#END OF DATA FORMATTING FOR PASS ON TO chi-sq-VR
results = pd.DataFrame(data = corr_final_2, columns = corr_bin_fin[1], index = list(range(1,(len(corr_final_2)+1))))
results_x = results.transpose()
results_x.to_csv("gc_correction_results", encoding='utf-8', index=True)
