# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 11:58:34 2025

@author: benco
"""

text = """
ADP                 1 GPVI = 1 ADP                         1
aggregation         5 Inta2b3 = 1 irrev.agg                3
Akt_activation      1 PI3 = 1 Akt                          0
Akt_activation2     1 Src + 2 PI3 = 2 Akt                  1
Akt_activation3     2 Src + 2 PI3 = 3 Akt                  3
ARAC_recruitment    1 PLA2 + 4 Inta2b3 = 1 ARAC            3
CaC                 1 IP3R = 1 CaC                         0
CaC2                2 IP3R + 1 P2X1 = 2 CaC                1
CaC3                1 Orai1 = 3 CaC                        2
CaC4                3 IP3R + 1 P2X1 + 1 Orai1 = 4 CaC      3
CalDAG_activation   1 CaC = 1 CalDAG                       0
CalDAG_activation2  2 CaC = 2 CalDAG                       1
CalDAG_activation3  3 CaC = 3 CalDAG                       2
CalDAG_activation4  4 CaC = 4 CalDAG                       3
cAMP_production     1 !P2Y12 + 1 PTGDR + 1 !PDE3 = 1 cAMP  0
cellshape           1 CaC + 1 RAP1 + 1 !VASP = 1 cellshape 0
cellshape2          2 CaC + 2 RAP1 = 1 cellshape           1
cellshape3          2 CaC + 2 RAP1 + 1 RhoK = 2 cellshape  3
DAG_production      1 PLCB = 1 DAG                         0
DAG_production2     2 PLCB = 2 DAG                         1
DAG_production3     3 PLCB = 3 DAG                         3
GPVI_activation      = 1 GPVI                              0
Int_activation      1 Tallin = 1 Inta2b3                   0
Int_activation1     1 Akt = 1 Inta2b3                      0
Int_activation1_2   1 Akt + 1 Tallin = 2 Inta2b3           0
Int_activation2     2 Akt + 2 Tallin = 3 Inta2b3           1
Int_activation3     2 Akt + 3 Tallin = 4 Inta2b3           2
Int_activation4     3 Akt + 3 Tallin = 5 Inta2b3           3
IP3_production      1 PLCB = 1 IP3                         0
IP3_production2     2 PLCB = 2 IP3                         1
IP3_production3     3 PLCB = 3 IP3                         3
IP3R_activation     1 PKC + 1 IP3 + 1 !PKA = 1 IP3R        0
IP3R_activation2    2 PKC + 2 IP3 = 2 IP3R                 1
IP3R_activation3    3 PKC + 3 IP3 = 3 IP3R                 3
Orai1_activation    2 CaC = 1 Orai1                        2
P2X1                1 ADP = 1 P2X1                         0
P2Y1                1 ADP = 1 P2Y1                         1
P2Y12_activation    1 ADP = 1 P2Y12                        1
PDE3_inhibition     1 PKA = 1 PDE3                         2
PI3_activation      1 IP3 = 1 PI3                          0
PI3_activation2     1 IP3 + 1 P2Y12 = 2 PI3                1
PKA_activation      1 cAMP = 1 PKA                         0
PKC_activation      1 DAG = 1 PKC                          0
PKC_activation2     1 CaC + 2 DAG = 2 PKC                  1
PKC_activation3     3 CaC + 2 DAG = 3 PKC                  2
PLA2_activation     1 P2Y1 = 1 PLA2                        1
PLCB_activation2    1 P2Y1 = 2 PLCB                        1
PLCB_activation3    1 TBXA2R = 3 PLCB                      3
PLCB_GPVI           1 GPVI = 1 PLCB                        0
PTGDR_activation     = 1 PTGDR                             0
RAP1_activation     1 CalDAG + 1 !PKA = 1 RAP1             0
RAP1_activation2    2 CalDAG + 1 PKA = 1 RAP1              1
RAP1_activation3    2 CalDAG + 1 !PKA = 2 RAP1             1
RAP1_activation4    3 CalDAG + 1 !PKA = 4 RAP1             2
RAP1_activation5    3 CalDAG + 1 PKA = 2 RAP1              2
RAP1_activation6    4 CalDAG = 4 RAP1                      3
RhoK_activation     1 TBXA2R + 1 !PKA = 1 RhoK             3
Src_activation      1 PKC + 1 !PKA = 1 Src                 0
Src_activation2     2 PKC = 1 Src                          1
Src_activation3     3 PKC = 2 Src                          2
Tallin_activation   1 RAP1 = 1 Tallin                      0
Tallin_activation2  2 RAP1 = 2 Tallin                      1
Tallin_activation3  3 RAP1 = 3 Tallin                      2
Tallin_activation4  4 RAP1 = 4 Tallin                      2
TBXA2R_activation   1 TXA = 1 TBXA2R                       3 
TXA_production      1 ARAC = 1 TXA                         3
VASP_activation     1 PKA = 1 VASP                         0
""".split('\n')

pmid = '23463387'

while '' in text: text.remove('')
for i, line in enumerate(text):
    line = line.split(' ')
    while '' in line: line.remove('')
    line = ''.join(line[j] for j in range(1, len(line) - 1))
    if line[0] == '=':
        text[i] = ''
        continue
    line = line.split('=')
    prefix = line[1][1:] + '=' + line[1][:1] + '::'
    line = line[0].split('+')
    line = ''.join(line[j][1:] + '=' + line[j][:1] + ' AND ' for j in range(len(line)))[:-5].replace('!', 'NOT ')
    text[i] = prefix + line
while '' in text: text.remove('')

rdict = dict()
for line in text:
    try:
        rdict[line.split('::')[0]] = '(' + rdict[line.split('::')[0]] + ') OR (' + line.split('::')[1] + ')'
    except KeyError:
        rdict[line.split('::')[0]] = line.split('::')[1]
sort_order = sorted(rdict)

g = open(pmid + '.txt', 'w')
for key in sort_order:
    g.write(key + ' :\t' + rdict[key] + '\n')
g.close()