#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 13:00:29 2025

@author: benco
"""

text = """
Wg
 2 
KWg({Slp, 1},{Ci[act], 1}) = 1 
KWg({Slp, 1},{Ci[act], 2}) = 2 
KWg({Slp, 1},{Ci[act], 2},{Nkd, 2}) = 2 
Fz 1 KFz({Wg, [1,2]}) = 1 
Dsh
 1 
Slp 1 
Nkd
 2 
KDsh({Fz, 1}) = 1 
KSlp({Dsh, 1}) = 1 
KNkd(∅) = 1 
KNkd({Dsh, 1}) = 2 
En 1 KEn({Dsh, 1}) = 1 
KEn({Dsh, 1},{En, 1}) = 1 
KEn({Dsh, 1},{Nkd, [1,2]}) = 1 
KEn({Dsh, 1},{En, 1},{Nkd, [1,2]}) = 1 
Hh
 1 
Ptc
 2 
Pka 2 
KHh({En, 1}) = 1 
KPtc(∅) = 1 
KPtc({Ci[act], [1,2]}) = 2 
KPka({Ptc, 1}) = 2 
Ci 1 KCi(∅) = 1 
Ci[act]
 2 
Ci[rep]
 1 
KCi[act]({Ci, 1}) = 1 
KCi[act]({Ci, 1},{Dsh, 1},{Pka, [1,2]}) = 1 
KCi[act]({Ci, 1},{Dsh, 1}) = 2 
KCi[rep]({Ci, 1},{Pka, 2}) = 1
""".replace('[act]', '@@').replace('[rep]', '##').split('\n')
pmid = '18956339'

while '' in text: text.remove('')

stext = []
for line in text:
    if 'K' in line:
        stext.append(line.split('K')[1])

r = []
for line in stext:
    node = line.split(' ')
    while '' in node: node.remove('')
    node = node[0].split('(')[0] + '=' + node[len(node) - 1] + ' :\t'
    regulators = line.split('(')[1].split(')')[0].replace(' ', '').split('}')
    while '' in regulators: regulators.remove('')
    rightnodes = ''
    for i, regnode in enumerate(regulators):
        regnode = regnode.replace(',{', '').replace('{', '')
        if '[' in regnode:
            states = regnode.split('[')[1]
            states = states[:len(states) - 1]
            regnode = regnode.split('[')[0] + states
        if not regnode == '∅':
            rule = ''
            rsplit = regnode.split(',')
            term_prefix = rsplit[0] + '='
            for j in range(1, len(rsplit)):
                rule += term_prefix + rsplit[j] + ' OR '
            rightnodes += '(' + rule[:-4] + ') AND ' # remove trailing ' OR '
        else:
            rightnodes += '12345' # add padding equivalent in length to a trailing ' AND '
    
    # combine the prefix and rule, remove the trailing ' AND ', and correct the placeholder titles
    r.append((node + rightnodes[:-5]).replace('@@', '[act]').replace('##', '[rep]'))

# combine functions that result in the same value, as well as adjusting for basal values
rdict = dict()
for rule in r:
    rsplit = rule.replace('Nkd=1', 'Nkd=0').replace('Nkd=2', 'Nkd=1').replace('Ptc=1', 'Ptc=0').replace('Ptc=2', 'Ptc=1').split(':\t')
    new_rule = rsplit[1]
    if 'Ci=1' in new_rule:
        new_rule = new_rule.replace('Ci=1', '')
    if new_rule.replace('(', '').replace(')', '') == '':
        continue
    new_rule = new_rule.replace('() AND ', '')
    try:
        rdict[rsplit[0]] = '(' + rdict[rsplit[0]] + ') OR (' + new_rule + ')'
    except KeyError:
        rdict[rsplit[0]] = new_rule
sort_order = sorted(rdict)

g = open(pmid + '.txt', 'w')
for key in sort_order:
    g.write(key + ':\t' + rdict[key] + '\n')
g.close()