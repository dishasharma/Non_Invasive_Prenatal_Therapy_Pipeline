#chrfractions.SeparatedStrands <- function(nipt_sample){
#  sapply(X = nipt_sample[[autosomal_chromosome_reads]], FUN = function(x) (rowSums(x) / sum(x)) / 2)
#}
##'@export
#chrfractions.CombinedStrands <- function(nipt_sample){
#  sapply(X = nipt_sample[[autosomal_chromosome_reads]], FUN = function(x) (rowSums(x) / sum(x)))
#}

def chrfractions(count_array):
    x_total = sum(count_array)
    x_chrfrac = []
    for i in count_array:
        x_curfrac = i/x_total
        x_chrfrac.append(x_curfrac)
    return x_chrfrac

import pandas as pd
import numpy as np
import re

x = pd.read_csv('gc_corr_chi_sq_results.csv')
x = x.transpose()
x_lis = x.values.tolist()
x = pd.DataFrame(data = x_lis[1:], columns = x_lis[0])
print(x)
read_list = x.values.tolist()
x_labels = x.columns.tolist()
#print(pd.DataFrame(data = x_labels))

x_l = list(range(1,23))
chrwise_l = []
xy_l = []
for xb in x_l:
    rex_chr = "^chr"+str(xb)+":.*$"
    rex_v = re.compile(rex_chr)
    temp = []
    for b in x_labels:
        if rex_v.match(b):
            temp.append(b)
    chrwise_l.append(temp)

#print(pd.DataFrame(data = chrwise_l))
sx_l = []
sx_l.append('X')
sx_l.append('Y')
chr_read_l = []

for xa in chrwise_l:
    tempr = []
    for i in xa:
        tempr.append(list(x[i]))
    chr_read_l.append(tempr)

chr_xy_l = []
for xb in sx_l:
    rex_chr = "^chr"+str(xb)+":.*$"
    rex_v = re.compile(rex_chr)
    temp = []
    for b in x_labels:
        if rex_v.match(b):
            temp.append(b)
    xy_l.append(temp)

#print(pd.DataFrame(data = xy_l))
for xa in xy_l:
    tempr = []
    for i in xa:
        tempr.append(list(x[i]))
    chr_xy_l.append(tempr)


#print(pd.DataFrame(data = chr_xy_l[0]))
chr_tran = []
for d in chr_read_l:
    chr_tran.append(np.transpose(d))
#print(pd.DataFrame(data = chr_tran[0]))

chr_xy_tran = []
for d in chr_xy_l:
    chr_xy_tran.append(np.transpose(d))
#print(pd.DataFrame(data = chr_xy_tran[0]))
fin_tabl = []

for tabl in chr_tran:
    temps = []
    for lit in tabl:
        sumr = sum(lit)
        temps.append(sumr)
    fin_tabl.append(temps)

xy_tabl = []

for tabl in chr_xy_tran:
    temps = []
    for lit in tabl:
        sumr = sum(lit)
        temps.append(sumr)
    xy_tabl.append(temps)

print(pd.DataFrame(data = xy_tabl))
xy_control = np.transpose(xy_tabl)
print(pd.DataFrame(data = xy_control))
condf = pd.DataFrame(data = fin_tabl)
control_df = condf.transpose()
control_df.to_csv("gc_corr_chi_sq_chrwise_results.csv", encoding='utf-8', index=True)
control_table = control_df.values.tolist()

print(control_df)
'''
sample_df = pd.DataFrame(pd.read_csv("sample_file.csv"))

x = sample_df.values.tolist()
sample_table = []
for i in x:
    sample_table.append(i[1:])
'''
#matchQC score calculation single sample and control
'''
for i in range(0,22):
    if(i!=13 and i!=18 and i!=21): 
        score = score + ((obfrac_s["*"][i] - obfrac_c["@"][i])**2)
'''
#matchQC score calculation single sample and multiple control
'''
table_sample_count = len(sample_table)
obfrac_s = []

for i in range(0,table_sample_count):
    y = []
    y = sample_table[i]
    obfrac_s.append(chrfractions(y))
'''

table_control_count = len(control_table)
obfrac_c = []

for i in range(0,table_control_count):
    x = control_table[i]
    obfrac_c.append(chrfractions(x))

#print(pd.DataFrame(data = obfrac_c))

match_qc_arr = []

for j in range(0,table_control_count):
    t_score = 0
    for k in range(0,table_control_count):
        for i in range(0,21):
            if(i!=12 and i!=17 and i!=20): 
                t_score = t_score + ((obfrac_c[j][i] - obfrac_c[k][i])**2)
    match_qc_arr.append(t_score)

print(pd.DataFrame(data = match_qc_arr))
matchqc_mean = np.mean(match_qc_arr)
matchqc_sd = np.std(match_qc_arr)
print(matchqc_mean)
print(matchqc_sd)
mqc_z = []
for i in match_qc_arr:
    mqc_z.append(((i-matchqc_mean)/matchqc_sd))
print(pd.DataFrame(data = mqc_z))
valid_controls = []
valid_reads = []
for ix in range(len(mqc_z)):
    if abs(mqc_z[ix]) < 2:
        ful_t = [*control_table[ix],*xy_control[ix]]
        valid_controls.append(ful_t)
        valid_reads.append(read_list[ix])

finvaldf = pd.DataFrame(data = valid_controls)
findf = finvaldf.transpose()
findf.to_csv("alldone_chrwise_results.csv", encoding='utf-8', index=True)
finreaddf = pd.DataFrame(data = valid_reads, columns = x_labels)
finrdf = finreaddf.transpose()
finrdf.to_csv("alldone_binwise_results.csv", encoding='utf-8', index=True)