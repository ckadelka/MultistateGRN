#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 09:58:27 2025

@author: benco
"""

text = """
EMF1 FT
0 1
1 0
AP1 FUL TFL1 EMF1 LFY
X X 0 0 2
0 0 0 1 1
0 0 1,2 0 1
0 0 1,2 1 0
0 1 1,2 0 1
0 1 0 1 1
0 2 0,1 0 2
0 2 0 1 1
0 2 2 0 1
0 2 1,2 1 1
0 1 1,2 1 0
1 X 0,1 0 2
1 X 2 1 0
1 X 0,1 1 1
1 X 2 0 1
2 X 2 1 1
2 X 0,1 X 2
2 X X 0 2
AP1 LFY AP2 EMF1 TFL1
X X X 0 0
2 X X 1 0
0,1 2 X 1 0
1 0,1 1 1 0
1 0,1 0 1 1
0 0,1 X 1 2
LFY EMF1
0 1
1,2 0
TFL1 SEP
0 1
1,2 0
AG LFY TFL1 FT AP1
2 X X X 0
0 X X 2 2
1 X X 2 1
0 LFY≥TFL1 0,1 2
1 LFY≥TFL1 0,1 1
0,1 LFY<TFL1 0,1 0
TFL1 AP2
0 1
1,2 0
AP1 TFL1 FUL
X 1,2 0
0 0 2
1 0 1
2 0 0
SEP LFY UFO AP3 PI AG AP1 AP3
1 X X 1,2 1,2 1,2 X 2
1 X X 1,2 1,2 X 1,2 2
0 1,2 1 X X X X 1
X 1,2 1 0 X X X 1
X 1,2 1 X 0 X X 1
X 1,2 1 X X 0 0 1
0 X 0 X X X X 0
X X 0 0 X X X 0
X X 0 X 0 X X 0
X X 0 X X 0 0 0
0 0 X X X X X 0
X 0 X 0 X X X 0
X 0 X X 0 X X 0
X 0 X X X 0 0 0
LFY SEP AP3 PI AG AP1 PI
X X 0 X 0 X 0
X 1 1,2 1,2 1,2 X 2
X 1 1,2 1,2 X 1,2 2
1,2 0 1,2 1,2 1,2 X 1
1,2 1 0 1,2 1,2 X 1
1,2 1 1,2 0 1,2 X 1
1,2 1 1,2 1,2 0 0 1
0 0 1,2 1,2 1,2 X 0
0 1 0 1,2 1,2 X 0
0 1 1,2 0 1,2 X 0
0 1 1,2 1,2 0 0 0
TFL1 LFY AP1 AG WUS AP2 SEP LUG CLF AG
TFL1≤LFY X X X 0 X X X 2
TLF1>LFY X X X X X X X 0
TLF1<LFY X X 1 1 X X X 2
2 2 X X 1 1 X X X 2
TFL1=LFY<2 X X X 1 X X X 0
TFL1<LFY X 0,1 0 1 X X 0 1
TFL1<LFY X 2 0 1 0 X 0 1
TFL1<LFY X 0,1 0 1 X 0 X 1
TFL1<LFY X 2 0 1 0 0 X 1
TFL1<LFY 0 0,1 0 1 X X X 1
TFL1<LFY 0 2 0 1 0 X X 1
2 2 X 0,1 0 1 X X 0 1
2 2 X 2 0 1 0 X 0 1
2 2 X 0,1 0 1 X 0 X 1
2 2 X 2 0 1 0 0 X 1
2 2 0 0,1 0 1 X X X 1
2 2 0 2 0 1 0 X X 1
TFL1<LFY 1,2 0,1 0 1 X 1 1 0
TFL1<LFY 1,2 2 0 1 0 1 1 0
2 2 1,2 0,1 0 1 X 1 1 0
2 2 1,2 2 0 1 0 1 1 0
TFL1<LFY X 2 0 1 1 X X 2
2 2 X 2 0 1 1 X X 2
WUS AG SEP WUS
0 X X 0
1 2 1 0
1 2 0 1
1 0,1 X 1
""".replace('≥', '>=').replace('≤', '<=').split('\n')
pmid = '15486106'

header_lines = []
for i, line in enumerate(text):
    if line == '':
        continue
    if not ('0' in line or 'X' in line) and not ('<' in line or '>' in line):
        header_lines.append(i)
header_lines.append(len(text))
r = []
for i in range(0, len(header_lines) - 1):
    if text[header_lines[i]] == '':
        continue
    nodes = text[header_lines[i]].split(' ')
    while '' in nodes: nodes.remove('')
    prefix = nodes[len(nodes) - 1] + '='
    j = header_lines[i] + 1
    j_esc = header_lines[i + 1]
    while j < j_esc:
        states = text[j].split(' ')
        while '' in states: states.remove('')
        if len(states) == 0:
            j += 1
            continue
        rule = prefix + states[len(states) - 1] + ' :/t '
        for k in range(len(states) - 1):
            if '>=' in states[k]:
                split = states[k].split('>=')
                states[k] = 'G'
                states.insert(k, 'G')
            elif '<=' in states[k]:
                split = states[k].split('<=')
                states[k] = 'L'
                states.insert(k, 'L')
            elif '=' in states[k]:
                split = states[k].split('=')
                states[k] = 'E'
                states.insert(k, 'E')
            elif '>' in states[k]:
                split = states[k].split('>')
                states[k] = 'g'
                states.insert(k, 'g')
            elif '<' in states[k]:
                split = states[k].split('<')
                states[k] = 'l'
                states.insert(k, 'l')
        for k in range(len(states) - 1):
            # X indicates it can be any state, and therefore means that
            # it is irelevant for this logical statement
            if not 'X' in states[k]:
                rule += nodes[k] + '=' + states[k] + ' AND '
        r.append(rule[:-5]) # remove the trailing AND
        j += 1

for i in range(len(r)):
    if '=G' in r[i]:
        nodes = []
        locations = []
        rsplit = r[i].split(' ')
        for j, x in enumerate(rsplit):
            if '=G' in x:
                nodes.append(x.split('=G')[0])
                locations.append(j)
        for _ in range(3):
            rsplit.pop(locations[0])
        # the only instances this occurs are for trinary nodes A >= B
        replacement = '((' + nodes[1] + '=0 AND (' + nodes[0] + '=0 OR ' + nodes[0] + '=1 OR ' + nodes[0] + '=2)) OR (' + nodes[1] + '=1 AND (' + nodes[0] + '=1 OR ' + nodes[0] + '=2)) OR (' + nodes[1] + '=2 AND ' + nodes[0] + '=2))' 
        rsplit.insert(locations[0], replacement)
        replacement = ''
        for _, v in enumerate(rsplit):
            replacement += v + ' '
        r[i] = replacement
    if '=g' in r[i]:
        nodes = []
        locations = []
        rsplit = r[i].split(' ')
        for j, x in enumerate(rsplit):
            if '=g' in x:
                nodes.append(x.split('=g')[0])
                locations.append(j)
        for _ in range(3):
            rsplit.pop(locations[0])
        # the only instances this occurs are for trinary nodes A > B
        replacement = '((' + nodes[1] + '=0 AND (' + nodes[0] + '=1 OR ' + nodes[0] + '=2)) OR (' + nodes[1] + '=1 AND ' + nodes[0] + '=2))'
        rsplit.insert(locations[0], replacement)
        replacement = ''
        for _, v in enumerate(rsplit):
            replacement += v + ' '
        r[i] = replacement
    if '=L' in r[i]:
        nodes = []
        locations = []
        rsplit = r[i].split(' ')
        for j, x in enumerate(rsplit):
            if '=L' in x:
                nodes.append(x.split('=L')[0])
                locations.append(j)
        for _ in range(3):
            rsplit.pop(locations[0])
        # the only instances this occurs are for trinary nodes A <= B
        replacement = '((' + nodes[0] + '=0 AND (' + nodes[1] + '=0 OR ' + nodes[1] + '=1 OR ' + nodes[1] + '=2)) OR (' + nodes[0] + '=1 AND (' + nodes[1] + '=1 OR ' + nodes[1] + '=2)) OR (' + nodes[0] + '=2 AND ' + nodes[1] + '=2))' 
        rsplit.insert(locations[0], replacement)
        replacement = ''
        for _, v in enumerate(rsplit):
            replacement += v + ' '
        r[i] = replacement
    if '=l' in r[i]:
        nodes = []
        locations = []
        rsplit = r[i].split(' ')
        for j, x in enumerate(rsplit):
            if '=l' in x:
                nodes.append(x.split('=l')[0])
                locations.append(j)
        for _ in range(3):
            rsplit.pop(locations[0])
        # the only instances this occurs are for trinary nodes A < B
        replacement = '((' + nodes[0] + '=0 AND (' + nodes[1] + '=1 OR ' + nodes[1] + '=2)) OR (' + nodes[0] + '=1 AND ' + nodes[1] + '=2))'
        rsplit.insert(locations[0], replacement)
        replacement = ''
        for _, v in enumerate(rsplit):
            replacement += v + ' '
        r[i] = replacement
    if '=E' in r[i]:
        nodes = []
        locations = []
        rsplit = r[i].split(' ')
        for j, x in enumerate(rsplit):
            if '=E' in x:
                nodes.append(x.split('=E')[0])
                locations.append(j)
        for _ in range(3):
            rsplit.pop(locations[0])
        # the only instances this occurs are for trinary nodes A=B<2
        replacement = '((' + nodes[0] + '=0 AND ' + nodes[1] + '=0) OR (' + nodes[0] + '=1 AND ' + nodes[1] + '=1))'
        rsplit.insert(locations[0], replacement)
        replacement = ''
        for _, v in enumerate(rsplit):
            replacement += v + ' '
        r[i] = replacement
    if ',' in r[i]:
        rsplit = r[i].split(' ')
        # the only instances this occurs are for a division between 2 states (i.e. 0,1; 0,2; 1,2)
        for j, x in enumerate(rsplit):
            if ',' in x:
                comma_split = x.split(',')
                replacement = '(' + comma_split[0] + ' OR ' +comma_split[0][0:len(comma_split[0]) - 1] + comma_split[1] + ')'
                rsplit[j] = replacement
        replacement = ''
        for _, v in enumerate(rsplit):
            replacement += v + ' '
        r[i] = replacement

# combine multiple lines that result in the same state for a given node
rdict = dict()
for _, rule in enumerate(r):
    rsplit = rule.split(':/t')
    try:
        rdict[rsplit[0]] = '(' + rdict[rsplit[0]] + ') OR (' + rsplit[1] + ')'
    except KeyError:
        rdict[rsplit[0]] = rsplit[1]
# remove excess spaces to achieve a nice look and slightly reduce file size
for key in rdict:
    rdict[key] = rdict[key].replace(' ', '').replace('AND', ' AND ').replace('OR', ' OR ')
sort_order = sorted(rdict)

g = open(pmid + '.txt', 'w')
for key in sort_order:
    g.write(key + ':\t' + rdict[key] + '\n')
g.close()