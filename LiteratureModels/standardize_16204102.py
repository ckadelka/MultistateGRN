#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  9 19:03:12 2025

@author: benco
"""

text = """
TrpE
TrpR = 0 and Trp = 0	1
TrpR = 1 and Trp = 0	0
TrpR = 0 and Trp ≥ 1	0
TrpR = 1 and Trp ≥ 1	0
TrpR
Trp ≤ 1	0
Trp = 2	1
Trp
Trpext = 0 and TrpE = 0	0
Trpext = 0 and TrpE = 1	1
Trpext = 1 and TrpE = 0	1
Trpext = 1 and TrpE = 1	1
Trpext = 2 and TrpE = 0	2
Trpext = 2 and TrpE = 1	2
""".replace(' ', '').split('\n')
pmid = '16204102'

r = []
for i, line in enumerate(text):
    if line == '':
        continue
    if not '\t' in line:
       current = line
    elif not line[len(line) - 1] == '0':
        prefix = current + '=' + line[len(line) - 1] + ' :\t'
        next_rule = line.split('\t')[0].replace('and', ' AND ')
        preexisting_rule = False
        for j, rule in enumerate(r):
            if prefix in rule:
                preexisting_rule = True
                if 'OR' in rule:
                    r[j] += ' OR (' + next_rule + ')'
                else:
                    r[j] = rule.split('\t')[0] + '\t(' + rule.split('\t')[1] + ') OR (' + next_rule + ')'
        if not preexisting_rule:
            r.append(prefix + next_rule)

g = open(pmid + '.txt', 'w')
for line in r:
    g.write(line + '\n')
g.close()