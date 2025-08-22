#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 20:52:19 2025

@author: benco
"""

text_b = """
TAB1 = !p38alpha
CFLAR = AKT & !ITCH
IKKA = AKT
DKK1 = DKK1gene
pras40 = !AKT
Axin = !LRP
TCF = betacatenin | !NLK
cMYC = TCF
RTPKgene = FOXO
TSC2 = GSK3 & !(IKKB | AKT | RSK | ERK)
SFRP1gene = !cMYC
Caspase9 = CytochromeC
p53 = p38alpha & !MDM2
DKK1gene = TCF & !cMYC
GAB = GRB2 & !ERK
AKT = PDK1 | mTORC2
ASK1 = !AKT
RSK = ERK & PDK1
SHP2 = GAB
Ras = SOS | SHP2
MEKK4 = Rac
S6K = PDK1 | mTORC1
MKK3 = ASK1 | TAK1
PDK1 = PI3K & !PTEN
MEK = Raf | MAP3K8 | !ERK
DUSP1 = p38alpha | MSK
BAD = !AKT & !RSK
BAX = p53
TAK1 = TAB1
RTPK = (RTPKgene | MMP) & !(p38alpha | MEK)
CK1 = !LRP
Egr1 = !TCF
SOS = GRB2 | !ERK
BCL2 = !BAD
MKK7 = TAK1 | GRAP2
LRP = (Fz | ERK | JNK1 | p38alpha) & !DKK1
GRB2 = SHC1
MAP3K8 = IKKB
Caspase8 = !CFLAR
FOXO = !(AKT | NLK)
GSK3 = !(LRP | RSK | S6K | ERK | p38alpha | Dvl | AKT)
Raf = Ras | !(Rheb | AKT | ERK)
ITCH = JNK1
MLK3 = Rac
PTENgene = Egr1
p38alpha = (MKK3 | MKK4) & !DUSP1
IKKB = TAK1 & !p53
MSK = ERK | p38alpha
MDM2 = (AKT | MDM2gene) & !S6K
MDM2gene = NFkB | p53
DUSP6 = ERK | mTORC1
NFkB = IKKA | IKKB | MSK
JNK1 = (MKK7 | MKK4) & !DUSP1
ERK = MEK | !DUSP6
Rheb = !TSC2
Rac = Dvl | mTORC2
CytochromeC = BAX & !BCL2
betacatenin = IKKA | !betaTrCP
MKK4 = MEKK4 | MLK3 | TAK1 | GRAP2
mTORC2 = TSC2 & !S6K
SHC1 = RTPK | !PTEN
IRS1 = !(S6K | ERK | IKKB)
mTORC1 = (Rheb | RSK) & !pras40
NLK = TAK1
Dvl = Fz
betaTrCP = Axin & GSK3 & CK1
SFRP1 = SFRP1gene
Fz = !SFRP1
MMP = LEF
PI3K = GAB | IRS1 | Ras
LEF = betacatenin
PTEN = PTENgene & !GSK3
GRAP2 = !p38alpha
""".split('\n')
text_m = """
CCND1 = RSK + TCF
Caspase37 = Caspase8 + Caspase9
Antisurvival = Caspase37 + FOXO
Prosurvival = CCND1 + cMYC
""".split('\n')
pmid = '26317215'

while '' in text_b: text_b.remove('')
while '' in text_m: text_m.remove('')

r = []
for i, line in enumerate(text_b):
    line = line.replace('!', ' NOT ').replace('|', ' OR ').replace('&', ' AND ').split('=')
    rule = line[0].replace(' ', '') + '=1 :\t'
    line = line[1].split(' ')
    while '' in line: line.remove('')
    for j, term in enumerate(line):
        if not (term == 'AND' or term == 'OR' or term == 'NOT'):
            if term[len(term) - 1] == ')':
                term = term[:len(term) - 1] + '=1)'
            else:
                term += '=1'
        rule += term + ' '
    r.append(rule[:len(rule) - 1])
for i, line in enumerate(text_m):
    line = line.split('=')
    node = line[0].replace(' ', '')
    line = line[1].replace(' ', '').split('+')
    while '' in line: line.remove('')
    # we know the rule for the multistate nodes consists of the sum of two other nodes
    # we also know that of the multistate nodes, the first two are ternary
    if i < 2:
        r.append(node + '=1 :\t(' + line[0] + '=1 AND ' + line[1] + '=0) OR (' + line[0] + '=0 AND ' + line[1] + '=1)')
        r.append(node + '=2 :\t' + line[0] + '=1 AND ' + line[1] + '=1')
    # and that the last two are quaternary, and furthermore that the quaternary nodes have a rule
    # consisting of a ternary node followed by a binary node
    else:
        r.append(node + '=1 :\t(' + line[0] + '=1 AND ' + line[1] + '=0) OR (' + line[0] + '=0 AND ' + line[1] + '=1)')
        r.append(node + '=2 :\t(' + line[0] + '=1 AND ' + line[1] + '=1) OR (' + line[0] + '=2 AND ' + line[1] + '=0)')
        r.append(node + '=3 :\t' + line[0] + '=2 AND ' + line[1] + '=1')

g = open(pmid + '.txt', 'w')
for line in r:
    g.write(line + '\n')
g.close()