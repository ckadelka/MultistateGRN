# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 13:22:48 2025

@author: benco
"""

text = """
S 1:
2:
ATM 1: S=1 AND E2F1=1 AND NOT (Wip1)
2: (S=2 OR E2F1=2) AND NOT (Wip1)
ATR 1: S=1
2: S=2
p53 1: (ATM=2 OR ATR=2 OR E2F1=2) AND (Mdm2 OR Wip1) AND NOT (Mdm2 AND Wip1)
2: (ATM OR ATR OR E2F1=2) AND NOT (Mdm2) AND NOT (Wip1)
miR16 1: p53=1 AND E2F1=1
2: p53=2 OR E2F1=2
Mdm2 1: (p53 OR pRB) AND NOT (ATM) AND NOT (Wip1)
Wip1 1: p53 AND NOT (miR16)
p21 1: p53 AND NOT (Mdm2)
pRB 1: NOT (Cdk2cE) AND NOT (Mdm2) AND NOT (Cdk46cD)
E2F1 1: (ATM OR ATR OR NOT (S)) AND NOT (p53) AND NOT (pRB) AND NOT (ATM=2 AND ATR=2)
2: ATM=2 AND ATR=2 AND NOT (p53) AND NOT (pRB)
Cdk2cE 1: NOT (p21) AND NOT (miR16) AND (E2F1 OR Cdc25)
Cdc25 1: NOT (ATM) OR NOT (ATR)
Cdk46cD 1: NOT (p21) AND NOT (miR16) AND Cdc25
""".split('\n')
pmid = '28968438'

while '' in text: text.remove('')

rdict = dict()
for i, line in enumerate(text):
    if line.split(':')[1] == '':
        continue
    node = line.split(' ')[0]
    if line[0] == '2':
        node = text[i - 1].split(' ')[0] + '=2'
    else:
        node += '=1'
    rule = line.split(':')[1][1:]
    rsplit = rule.split(' ')
    for j, r in enumerate(rsplit):
        if r == 'AND' or r == 'OR' or r == 'NOT' or '=' in r:
            continue
        if r.replace('(', '').replace(')', '') in [ 'Mdm2', 'Wip1', 'p21', 'pRB', 'Cdk2cE', 'Cdc25', 'Cdk46cD' ]:
            if r[len(r) - 1] == ')':
                rsplit[j] = r[:len(r) - 1] + '=1)'
            else:
                rsplit[j] = r + '=1'
    rdict[node] = ''.join(rsplit[j] + ' ' for j in range(len(rsplit))).replace('NOT (S)', 'S=0').replace('NOT (ATM)', 'ATM=0').replace('NOT (ATR)', 'ATR=0').replace('NOT (p53)', 'p53=0').replace('NOT (miR16)', 'miR16=0').replace('NOT (E2F1)', 'E2F1=0').replace('(ATM OR ATR OR E2F1=2)', '(ATM=2 OR ATR=2 OR E2F1=2)').replace('(p53 OR pRB=1)', '(NOT (p53=0) OR pRB=1)').replace('p53 AND', 'NOT (p53=0) AND').replace('(ATM OR ATR OR S=0)', '(NOT (ATM=0) OR NOT (ATR=0) OR S=0)').replace('(E2F1 OR Cdc25=1)', '(NOT (E2F1=0) OR Cdc25=1)')

g = open(pmid + '.txt', 'w')
for key in rdict.keys():
    g.write(key + ' :\t' + rdict[key] + '\n')
g.close()