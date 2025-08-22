#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  2 17:28:26 2025

@author: benco
"""

text = """
Sry 
1 
KSry
 ({Y,1},{Wt1,1},{AS,1},{Sf1,1},{Gata4,1})=1 
Sf1 2 KSf1
 ({Wt1,1})=1 
KSf1
 ({Wt1,1},{Sox9,2})=2 
Sox9 
2 
Fgf9 2 
KSox9
 ({Sf1,1})=1 
KSox9
 ({Sf1,1}, {b-cat,1})=1  
KSox9
 ({Sf1,1}, {b-cat,1},{Sox9,1})=1  
KSox9
 ({Sf1,1},{Sry,1})=2 
KSox9
 ({Sf1,1},{b-cat,1},{Sry,1})=2 
KSox9
 ({Sf1,1},{Sry,1},{Sox9,1})=2 
KSox9
 ({Sf1,1},{b-cat,1},{Sry,1},{Sox9,1})=2 
KSox9
 ({Sf1,1},{Sox9,1})=2 
KFgf9
 ({Fgf9,1})=1 
KFgf9
 ({Sox9,2})=1 
KFgf9
 ({Fgf9,1},({Sox9,2})=2 
Wnt4 2 KWnt4
 ({Wnt4,1},{Fgf9,1})=1 
KWnt4
 ({Wnt4,1})=2 
b-cat 
2 Kb-cat
 ({Wnt4,1})=1  
Kb-cat
 ({Wnt4,2})=2 
Dmrt1 
1 
KDmrt
 ({Gata4,1})=1  
KDmrt
 ({Gata4,1},{Sf1,1})=1  
KDmrt
 ({Gata4,1},{Sox9,2})=1  
KDmrt
 ({Gata4,1},{Sox9,2},{Sf1,1})=1  
KDmrt
 ({Sox9,2},{Sf1,1})=1  
Foxl2 
1 
Fgf9_c 2 
KFoxl2
 ({AF,1})=1 
KFgf9_c
 ({Sox9_c,2})=1 
KFgf9_c
 ({Fgf9r_c,1})=1 
KFgf9_c
 ({Sox9_c,2},{Fgf9r_c,1})=2 
Fgf9r_c 2 KFgf9r_c
 ({Fgf9_c,1})=1 
KFgf9r_c
 ({Fgf9_c,1},({Fgf9_p,2})=2 
Wnt4_c 2 KWnt4_c
 ({Wnt4_c,1},{Fgf9r_c,1})=1 
KWnt4_c
 ({Wnt4_c,1})=2  
Fgf9_p 2 
KFgf9_p
 ({Sox9_p,2})=1 
KFgf9_p
 ({Fgf9r_p,1})=1 
KFgf9_p
 ({Sox9_p,2},{Fgf9r_p,1})=2 
Fgf9r_p 2 KFgf9r_p
 ({Fgf9_p,1})=1 
KFgf9r_p
 ({Fgf9_p,1},({Fgf9_c,2})=2 
Wnt4_p 2 2KWnt4_p
 ({Wnt4_p,1},{Fgf9r_p,1})=1 
KWnt4_p
 ({Wnt4_p,1})=2
""".split('\n')
pmid = '27229461'

for i in range(len(text)):
    if 'K' in text[i]:
        text[i] = text[i].split('K')[1].replace(' ', '')
    elif '(' in text[i]:
        text[i - 1] += ' ' + text[i].replace(' ', '')
        text[i] = ''
    else:
        text[i] = ''
while '' in text: text.remove('')

r = []
for line in text:
    node = line.split(' ')
    while '' in node: node.remove('')
    regulators = node[1]
    node = node[0] + '=' + regulators.split('=')[1] + ' :\t'
    regulators = regulators.split('=')[0].replace(')', '').replace('(', '').split('}')
    while '' in regulators: regulators.remove('')
    rightnodes = ''
    for i, regnode in enumerate(regulators):
        regnode = regnode.replace(',{', '').replace('{', '')
        rsplit = regnode.split(',')
        rule = rsplit[0] + '=' + rsplit[1]
        rightnodes += rule + ' AND '
    r.append((node + rightnodes[:-5]))

rdict = dict()
for rule in r:
    node = rule.split(':\t')
    new_rule = node[1]
    node = node[0]
    try:
        rdict[node] = '(' + rdict[node] + ') OR (' + new_rule + ')'
    except KeyError:
        rdict[node] = new_rule
sort = sorted(rdict)

g = open(pmid + '.txt', 'w')
for key in sort:
    g.write(key + ':\t' + rdict[key] + '\n')
g.close()