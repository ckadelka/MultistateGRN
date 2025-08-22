#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 13:51:37 2025

@author: benco
"""

text = """
CtrAat	GcrAt	CcrMt	SciPt	CtrAa t+1
0	0	0	0	1
0	0	0	1	0
0	0	1	0	0
0	0	1	1	0
0	1	0	0	1
0	1	0	1	0
0	1	1	0	0
0	1	1	1	0
1	0	0	0	1
1	0	0	1	0
1	0	1	0	0
1	0	1	1	0
1	1	0	0	1
1	1	0	1	0
1	1	1	0	0
1	1	1	1	0
CtrAat	DnaAt	GcrA t+1
0	0	0
0	1	1
1	0	0
1	1	0
CtrAat	GcrAt	DnaAt	CcrMt	DnaA t+1
0	0	0	0	0
0	0	0	1	0
0	0	1	0	0
0	0	1	1	0
0	1	0	0	0
0	1	0	1	0
0	1	1	0	0
0	1	1	1	0
1	0	0	0	0
1	0	0	1	1
1	0	1	0	0
1	0	1	1	0
1	1	0	0	0
1	1	0	1	0
1	1	1	0	0
1	1	1	1	0
CtrAat	CcrMt	SciPt	CcrM t+1
0	0	0	0
0	0	1	0
0	1	0	0
0	1	1	0
1	0	0	1
1	0	1	0
1	1	0	0
1	1	1	0
DnaAt	CtrAat	SciP t+1
0	0	0
0	1	1
1	0	0
1	1	0
ChpTt	ClpXP-RcdAt	CtrAbt+1
1	0	1
1	1	0
2	0	2
2	1	0
PleCt	DivJt	DivK t+1
1	1	1
1	2	2
2	1	1
2	2	1
DivKt	PleC t+1
2	1
1	2
DivKt	PleCt	DivJ t+1
1	1	1
2	1	2
1	2	1
2	2	1
DivKt	DivL t+1
1	1
2	0
DivLt	CckA t+1
1	2
0	1
CckAt	ChpT t+1
2	2
1	1
ChpTt	CpdR t+1
2	2
1	1
CpdRt	ClpXP-RcdA t+1
1	1
2	0
""".split('\n')
pmid = "25369202"

r = []
for i, line in enumerate(text):
    if line == '':
        continue
    if 't+1' in line:
        line2 = line.split('t+1')[0].replace('t\t', ' ')
        nodes = line2.split(' ')
        while '' in nodes: nodes.remove('')
        prefix = nodes[len(nodes) - 1] + '='
        j = i + 1
        while j < len(text) - 1 and not 't+1' in text[j]:
            states = text[j].split('\t')
            while '' in states: states.remove('')
            rule = prefix + states[len(states) - 1] + ' :\t'
            for k in range(len(states) - 1):
                rule += nodes[k] + '=' + states[k] + ' AND '
            r.append(rule[:-5])
            j += 1
rdict = dict()
for rule in r:
    rsplit = rule.split(':\t')
    try:
        if 'OR' in rdict[rsplit[0]]:
            rdict[rsplit[0]] = rdict[rsplit[0]] + ' OR (' + rsplit[1] + ')'
        else:
            rdict[rsplit[0]] = '(' + rdict[rsplit[0]] + ') OR (' + rsplit[1] + ')'
    except KeyError:
        rdict[rsplit[0]] = rsplit[1]
sort_order = sorted(rdict)

g = open(pmid + '.txt', 'w')
for key in sort_order:
    g.write(key + ':\t' + rdict[key] + '\n')
g.close()