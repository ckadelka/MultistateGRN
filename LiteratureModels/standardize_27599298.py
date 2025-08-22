# -*- coding: utf-8 -*-
"""
Created on Wed Apr  2 17:46:02 2025

@author: benco
"""

text = """
 Bap 1 (Bin:1 | (Mad & Med)) & Tin & !Slp & !(Bin:1 & Mad & Med & Tin & En) & !(Ci & Tin & En) & !(Ci & Mad & Med & Tin & Bin:1) & !Bin:2
 1 Slp:1 & Mad:1 & Med & Ci & En & Tin & !Bin:2
 1 Ci & Tin & En & Slp & !(Mad & Med) & !Bin:2
 1 Ci & Tin & !En & Slp & Mad & Med & !Bin:2
 2 Bin:1 & Mad & Med & !(Ci & En) & !Slp & Tin & (Ci | En)
 2 Ci & Tin & En & !Slp & !(Mad & Med) & !Bin:2
 3 Bin:2 | (Ci:1 & En:1 & !Slp & Mad & Med & Tin)
 Bin 1 Bap
 Brk 1 !(Shn & Mad & Med)
 Ci 1 !Pka
 Da 1 !Emc
 Dsix4 1 (Tin:1 | Zfh-1:1) & !(Mad & Med)
 Dome 1 Upd
 Der 1 Spi:1
 Doc 1 Mad:1 & Med & Pan:1
 Emc 1 Nicd
 E_Spl 1 Nicd & Stat92E & SuH
 Eve 1 Pnt:1 & Htl:1 & Mad:1 & Pan:1 & Med & Tin:2 & Twi
 Eya 1 Twi | Tin:2
 Hbr 1 Pan & Htl
 Htl 1 (((Pyr & Pan) | Ths) & Pan) | (Pyr & Ths)
 Hop 1 Dome
 Mad 1 Tkv
 Mef2 1 Tin | Twi:2
 Nicd 1 Notch
 Notch 1 Delta
 Pan 1 Wg:1
 Pka 1 !Smo
 Pnr 1 Doc:1 & Tin
 Pnt 1 Rl:1
 Poxm 1 Pan & Twi & !(Mad & Med)
 Ptc 1 !Hh
 Ras 1 Der:1
 Rl 1 Ras:1
 Slp 1 Pan:1
 Smo 1 !Ptc
 Srp 1 Ci & !En & !(Mad & Med)
 2 Ci & En & !(Mad & Med)
 Stat92E 1 Hop | Tin:2
 Tin 1 Mad & (Med | Twi | Tin) & !(Mad & Med & Pan & Stat92E)
 2 Mad & Med & Pan & Stat92E
 Tkv 1 Dpp:1
 Twi 1 (Slp | Da) & Twi & !(Da & Slp) & !E_Spl
 1 Slp & Twi & !(Da & !E_Spl)
 2 Slp & Da & Twi & !E_Spl
 Zfh-1 1 Twi
""".split('\n')
pmid = '27599298'

for i, line in enumerate(text):
    if line == '':
        continue
    spline = line.split(' ')
    while '' in spline: spline.remove('')
    if spline[0] == '1' or spline[0] == '2' or spline[0] == '3':
        spline.insert(0, text[i - 1].split(' ')[0])
        spline[1] += '@@'
    elif spline[1] == '1' or spline[1] == '2' or spline[1] == '3':
        spline[1] += '@@'
    text[i] = ''.join(spline[j] + ' ' for j in range(len(spline)))
while '' in text: text.remove('')

r = []
for i, line in enumerate(text):
    node = line.split('@@')[0].replace(' ', '=')
    rule_terms = line.split('@@')[1]
    
    rule_terms = rule_terms.replace('(', ' ( ').replace(')', ' ) ').replace('!', ' NOT ').replace('&', ' AND ').replace('|', ' OR ').split(' ')
    while '' in rule_terms: rule_terms.remove('')
    
    for j, term in enumerate(rule_terms):
        if not (term == 'AND' or term == 'OR' or term == 'NOT' or term == '(' or term == ')'):
            if ':' in term:
                term = term.replace(':', '=')
            else:
                term += '=1'
        rule_terms[j] = term
    r.append(node + ' :\t' + ''.join(x + ' ' for x in rule_terms))
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
    g.write(key + ':\t' + rdict[key].replace('( ', '(').replace(' )', ')') + '\n')
g.close()