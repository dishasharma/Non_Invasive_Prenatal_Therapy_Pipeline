import os
import pandas as pd
import numpy as np
import glob
import math
import re

'''
maindir = 'bedfiles/'                                                  #functional directory here
dir_list = os.listdir(maindir)
#xdir = 'xyzabc.kmerlists'
#os.makedirs(xdir)
listing = glob.glob(os.path.join(maindir, '*.bed'))               #creates listing of bed files in directory
listing.sort()
final_read_list = []
final_bin_list = []

for fname in listing:
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
    final_read_list.append(content_readcount)

final_bin_list = content_chrname

x_df = pd.DataFrame(data = final_read_list, index = list(range(1,37)), columns = final_bin_list)
'''
rex_sample = "^chr[0-9]{1,2}.*$"
rex_v = re.compile(rex_sample)
x = pd.read_csv('gc_correction_results.csv')
csv_data_raw = x.values.tolist()
t_csv_dat = np.array(csv_data_raw).T.tolist()
final_bin_list = t_csv_dat[0]
final_read = t_csv_dat[1:]
sample_total_reads = []
total_reads = 0
final_read_list = []
auto_bin_list = []

for x in final_read:
    temp = []
    for y in x:
        tem_f = float(y)
        temp.append(tem_f)
    final_read_list.append(temp)

for x in final_bin_list:
    if rex_v.match(x):
        auto_bin_list.append(x)

bin_count = len(final_bin_list)
sample_count = len(final_read_list)

frdf = pd.DataFrame(data = final_read_list)
frdf = frdf.fillna(0)
final_read_list = frdf.values.tolist()

for x in final_read_list:
    sample_total_reads.append(np.nansum(x))

print(sample_total_reads)

total_reads = sum(sample_total_reads)

norm_factor_arr = []

for i in sample_total_reads:
    n_f = ((total_reads/(bin_count*sample_count))/(i/(bin_count)))
    norm_factor_arr.append(n_f)

norm_read_arr = []

for i in range(len(final_read_list)):
    temp = []
    te_st = final_read_list[i]
    for j in range(len(te_st)):
        x_val = te_st[j]*norm_factor_arr[i]
        temp.append(x_val)
    norm_read_arr.append(temp)

bin_wise_avg = []
transp_list = []
y_df = pd.DataFrame(data = norm_read_arr, columns = final_bin_list)

for wd in final_bin_list:
    hs = list(y_df[wd])
    transp_list.append(hs)
    loc_tot = np.nansum(hs)
    avrg = loc_tot/len(hs)
    bin_wise_avg.append(avrg)

#print(sample_total_reads)
#print(total_reads)
print(sample_count)
print(bin_count)
print(pd.DataFrame(data = norm_factor_arr,index = list(range(1,sample_count+1))))
#print(pd.DataFrame(data = bin_wise_avg,index = final_bin_list))
#for i in control_table:
chi_sq_sum = []
for x in range(len(transp_list)):
    sum = 0
    for u in transp_list[x]:
        sum = sum + (((bin_wise_avg[x] - u)**2)/bin_wise_avg[x])
    chi_sq_sum.append(sum)

norm_chi_sq = []
for i in chi_sq_sum:
    csd = (i - sample_count)/(math.sqrt(2)*sample_count)
    norm_chi_sq.append(csd)

df_norm = pd.DataFrame(data = norm_chi_sq, index = final_bin_list)
#df_norm = df_norm.fillna(0)
print(df_norm)
#to calculate normalization factor in chi-squared-based-VR
'''
pos = 0
neg = 0
zex = 0
nen = 0
for i in norm_chi_sq:
    if i > 0:
        pos = pos + 1
    elif i < 0:
        neg = neg + 1
    elif math.isnan(i):
        nen = nen + 1
    else:
        zex = zex + 1
print(pos)
print(neg)
print(zex)
print(nen)
'''
abs_chi_sq = []
corr_facc = []
for i in norm_chi_sq:
    if i > 0:
        abs_chi_sq.append(i)
    elif i < 0:
        abs_chi_sq.append(abs(i))
    elif math.isnan(i):
        abs_chi_sq.append(i)
    else:
        abs_chi_sq.append(i)

for j in range(len(abs_chi_sq)):
    if abs_chi_sq[j] > 2:
        corr_fac = chi_sq_sum[j]/sample_count
    else:
        corr_fac = 1
    corr_facc.append(corr_fac)

corr_read_list = []

print(pd.DataFrame(data = corr_facc, index = final_bin_list))
for s_l in final_read_list:
    temp = []
    for i in range(len(s_l)):
        corr_val = s_l[i]/corr_facc[i]
        temp.append(corr_val)
    corr_read_list.append(temp)

x = pd.DataFrame(data = corr_read_list, columns = final_bin_list)
y = pd.DataFrame(data = final_read_list, columns = final_bin_list)
#print(x['chr1:600000-650000'])
#print(y['chr1:600000-650000'])
results_x = x.transpose()
results_x.to_csv("gc_corr_chi_sq_results.csv", encoding='utf-8', index=True)