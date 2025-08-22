# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 14:56:04 2025

@author: benco
"""

text_core = """
Cln3 1 ● MASS
Bck2 1 ● MASS
MBF_SBF 1 ● !(Clb2 & !(Sic1 | Cdc6)) & (Bck2 | Cln3 | Cln2 | (Clb5 
& !Sic1))
Cln2 1 ● MASS & MBF_SBF
Swi5 1 ● (Mcm1 & !(Clb2 & !(Sic1 | Cdc6))) | (Mcm1 & Clb2 
& !(Sic1 | Cdc6) & ((Cdc14:1 & !Net1) | (Cdc14:2 & !
 Net1:2)))
Sic1 1 ● (!Cdc14 | (Cdc14:1 & Net1) | ((Cdc14:2 | Cdc14:3) & 
Net1:2)) & !Swi5 & !((Clb2 & !(Sic1 | Cdc6)) | (Clb5 
& !Sic1) | Cln2 | ((Clb2 | Clb5) & (Cln3 | Bck2:2)) | 
(Clb5 & Clb2) | (Clb5:3 & Bck2))
● (!Cdc14 | (Cdc14:1 & Net1) | ((Cdc14:2 | Cdc14:3) & 
Net1:2)) & Swi5 & !((Clb2:3 & !(Sic1 | Cdc6)) | 
(((Clb2 & !(Sic1 | Cdc6)) | (Clb5 & !Sic1)) & ((Cln2 & 
(Cln3 | Bck2)) | (Cln3 & Bck2))) | (((Clb2 & Clb5) | 
Clb2:3 | Clb5:3) & Cln2 & (Cln3 | Bck2)) | (Clb2:3 & 
Clb5:3))
● ((Cdc14:1 & !Net1) | (Cdc14:2 & (!Net1 | Net1:1))) 
& !Swi5 & !((Clb2 & !(Sic1 | Cdc6)) | (Clb5 & !Sic1) | 
Cln2)
● ((Cdc14:1 & !Net1) | (Cdc14:2 & (!Net1 | Net1:1))) & 
Swi5 & !(((Clb5 & Clb2 & !(Sic1 | Cdc6) & Cln2 & 
Cln3 & Bck2) | (Clb2:3 & !(Sic1 | Cdc6))) & (Clb5 | 
Cln2 | (Cln3 & Bck2)))
 ●Cdc14:3 & !Net1:2 & !Swi5 & !(Clb5 & Clb2 & Cln2 
& Cln3 & Bck2)
 ●Cdc14:3 & !Net1:2 & Swi5 & !(Cln3 & Bck2 & Cln2 
& Clb5 & !Sic1 & Clb2 & !(Sic1 | Cdc6))
 Cdc6 1 ● (!Cdc14 | (Cdc14:1 & Net1) | ((Cdc14:2 | Cdc14:3) & 
Net1:2)) & !Swi5 & !((Clb2 & !(Sic1 | Cdc6)) | (Clb5 
& !Sic1) | Cln2 | ((Clb2 | Clb5) & (Cln3 | Bck2:2)) | 
(Clb5 & Clb2) | (Clb5:3 & Bck2))
● (!Cdc14 | (Cdc14:1 & Net1) | ((Cdc14:2 | Cdc14:3) & 
Net1:2)) & Swi5 & !((Clb2:3 & !(Sic1 | Cdc6)) | 
(((Clb2 & !(Sic1 | Cdc6)) | (Clb5 & !Sic1)) & ((Cln2 & 
(Cln3 | Bck2)) | (Cln3 & Bck2))) | (((Clb2 & Clb5) | 
Clb2:3 | Clb5:3) & Cln2 & (Cln3 | Bck2)) | (Clb2:3 & 
Clb5:3))
 ● ((Cdc14:1 & !Net1) | (Cdc14:2 & (!Net1 | Net1:1))) 
& !Swi5 & !((Clb2 & !(Sic1 | Cdc6)) | (Clb5 & !Sic1) | 
Cln2)
● ((Cdc14:1 & !Net1) | (Cdc14:2 & (!Net1 | Net1:1))) & 
Swi5 & !(((Clb5 & Clb2 & !(Sic1 | Cdc6) & Cln2 & 
Cln3 & Bck2) | (Clb2:3 & !(Sic1 | Cdc6))) & (Clb5 | 
Cln2 | (Cln3 & Bck2)))
 ●Cdc14:3 & !Net1:2 & !Swi5 & !(Clb5 & Clb2 & Cln2 
& Cln3 & Bck2)
●Cdc14:3 & !Net1:2 & Swi5 & !(Cln3 & Bck2 & Cln2 
& Clb5 & !Sic1 & Clb2 & !(Sic1 | Cdc6))
 Clb5 1 ●MASS & MBF_SBF & Cdc20
2 ●MASS & MBF_SBF & !Cdc20
 Clb2 2 ●MASS & !Cdh1 & ((!Cdc20 & !Mcm1) | (Cdc20:2 & 
Mcm1))
 3 ●MASS & !Cdh1 & !Cdc20 & Mcm1
 Mcm1 1 ●Clb2 & !(Sic1 | Cdc6)
Mad2 1 ●ORI & !SPN
Cdc20 1 ● !Mad2 & Mcm1 & !(Clb2 & !(Sic1 | Cdc6))
 2 ● !Mad2 & Mcm1 & Clb2 & !(Sic1 | Cdc6)
PPX 1 ●Pds1 & !Cdc20
 Bub2Bfa1 1 ●ORI & !SPN
 Lte1 1 ● (Clb2 & !(Sic1 | Cdc6)) | SPN
Tem1 1 ●Lte1
Cdc15 1 ● (basal value)
Net1 1 ● !(((Cdc15:1 & Tem1) | Cdc15:2) & !Bub2Bfa1) | PPX
 Cdc14 1 ● (basal value)
 Cdh1 1 ● (!Cdc14 | (Cdc14:1 & Net1) | ((Cdc14:2 | Cdc14:3) & 
Net1:2)) & !((Clb2 & !(Sic1 | Cdc6)) | (Clb5 & !Sic1) | 
(Cln3 & Cln2))
 ● ((Cdc14:1 & !Net1) | (Cdc14:2 & !Net1:2)) & !((Clb5 
& !Sic1 & Cln3 & ((Clb2 & !Cdc6) | Cln2)) | (Clb2:3 
& !(Sic1 | Cdc6) & Cln3 & Cln2))
●Cdc14:3 & !Net1:2 & !((Clb5 & !Sic1 & Cln3 & Clb2 
& !Cdc6 & Cln2) | (Clb2:3 & !(Sic1 | Cdc6) & Cln3 & 
Cln2))
 BUD 1 ● (Cln2 | (Clb5 & !Sic1)) & !CYTOKINESIS
 ORI 1 ● (Clb5 & !Sic1) | ((Clb2:2 | Clb2:3) & !(Sic1 | Cdc6)) | 
(ORI & (Clb5:3 | Clb2:3 | (Clb5:1 & !Sic1) | (Clb2:1 & 
!(Sic1 | Cdc6))))
 SPN 1 ●Clb2 & !(Sic1 | Cdc6) & !CYTOKINESIS
 Pds1 1 ● (Mcm1 & MBF_SBF & !Cdc20) | ((Mcm1 | 
MBF_SBF) & !Cdh1 & !Cdc20)
 Esp1 1 ● !Pds1
 MASS 1 ● !CYTOKINESIS
 CYTOKINESIS 1 ● MASS & Clb2:2 & !(Sic1 | Cdc6)
 2 ● MASS & ((Clb2:1 & CYTOKINESIS) | (Clb2:2 & 
CYTOKINESIS & (Sic1 | Cdc6)))
""".split('\n')
while '' in text_core: text_core.remove('')

text_morpho = """
 BUD 1 ● MASS
 SBF 1 ● MASS & !Clb2
 Swe1 1 ●SBF & ((Clb2 & !Hsl1) | (Hsl1 & !
Clb2))
 2 ●SBF & !(Hsl1 | Clb2)
 Mih1 1 ● (Mpk1 & Clb2) | (!Mpk1 & !Clb2)
 2 ● !Mpk1 & Clb2
 Clb2 1 ● (MASS:1 & ((Swe1:1 & !Mih1:2) | 
(Swe1:2 & Mih1:1))) | (MASS & 
Swe1:2 & !Mih1)
2 ● (MASS:1 & (!Swe1 | Mih1:2)) | 
(MASS:2 & (!Swe1:2 | Mih1))
 MASS 1 ●MASS:1
 2 ●MASS:2
 Mpk1 1 ● !BUD
 Hsl1 1 ●BUD
""".split('\n')
while '' in text_morpho: text_morpho.remove('')

text_exit = """
Clb2 1 ● !Cdh1 & Cdc20
2 ● !Cdc20 & !Cdh1
Cdc20 1 ● !Cdh1
SecurinPds1 1 ● !Cdc20
SeparaseEsp1 1 ● (basal value)
PP2ACdc55 1 ● !(!SeparaseEsp1 | (SeparaseEsp1:1 & 
SecurinPds1))
2 ● !SeparaseEsp1 | (SeparaseEsp1:1 & 
SecurinPds1)
 Cdc5Polo 1 ●Clb2 & !Cdh1
 Bub2Bfa1 1 ●PP2ACdc55 | !Cdc5Polo
 Tem1 1 ● (basal value)
 Cdc15 1 ● !Clb2 | (Cdc14 & !Net1)
 Net1 1 ● ((Cdc14 & Net1 & !Clb2 & !
PP2ACdc55) | (PP2ACdc55:1 & ((!
Cdc14 & Clb2:1) | (Cdc14 & Clb2:1)))) & !(Cdc15 & Tem1 & !
 Bub2Bfa1)
2 ● ((((Cdc14 & !Net1) | PP2ACdc55) 
& !Clb2) | PP2ACdc55:2) & !(Cdc15 
& Tem1 & !Bub2Bfa1)
Cdc14 1 ● (basal value)
 Cdh1 1 ● (!Clb2:2 & Cdc14 & !Net1) | (Cdh1 
& !Clb2)
""".split('\n')
while '' in text_exit: text_exit.remove('')

text_coupled = """
Cln3 1 ● MASS
Bck2 1 ● MASS
MBF_SBF 1 ● !(Clb2 & !(Sic1 | Cdc6)) & (Bck2 | Cln3 | Cln2 | 
(Clb5 & !Sic1))
Cln2 1 ● MASS & MBF_SBF
 Swi5 1 ● (Mcm1 & !(Clb2 & !(Sic1 | Cdc6))) | (Mcm1 & 
Clb2 & !(Sic1 | Cdc6) & ((Cdc14:1 & !Net1) | 
(Cdc14:2 & !Net1:3)))
 Sic1 1 ● (!Cdc14 | (Cdc14:1 & Net1) | ((Cdc14:2 | 
Cdc14:3) & Net1:3)) & !Swi5 & !((Clb2 & !
 (Sic1 | Cdc6)) | (Clb5 & !Sic1) | Cln2 | ((Clb2 | 
Clb5) & (Cln3 | Bck2:2)) | (Clb5 & Clb2) | 
(Clb5:3 & Bck2))
● (!Cdc14 | (Cdc14:1 & Net1) | ((Cdc14:2 | 
Cdc14:3) & Net1:3)) & Swi5 & !((Clb2:3 & !
 (Sic1 | Cdc6)) | (((Clb2 & !(Sic1 | Cdc6)) | (Clb5 
& !Sic1)) & ((Cln2 & (Cln3 | Bck2)) | (Cln3 & 
Bck2))) | (((Clb2 & Clb5) | Clb2:3 | Clb5:3) & 
Cln2 & (Cln3 | Bck2)) | (Clb2:3 & Clb5:3))
 ● ((Cdc14:1 & !Net1) | (Cdc14:2 & (!Net1 | 
Net1:1))) & !Swi5 & !((Clb2 & !(Sic1 | Cdc6)) | 
(Clb5 & !Sic1) | Cln2)
 ● ((Cdc14:1 & !Net1) | (Cdc14:2 & (!Net1 | 
Net1:1))) & Swi5 & !(((Clb5 & Clb2 & !(Sic1 | 
Cdc6) & Cln2 & Cln3 & Bck2) | (Clb2:3 & !
 (Sic1 | Cdc6))) & (Clb5 | Cln2 | (Cln3 & Bck2)))
 ●Cdc14:3 & !Net1:3 & !Swi5 & !(Clb5 & Clb2 & 
Cln2 & Cln3 & Bck2)
 ●Cdc14:3 & !Net1:3 & Swi5 & !(Cln3 & Bck2 & 
Cln2 & Clb5 & !Sic1 & Clb2 & !(Sic1 | Cdc6))
 Cdc6 1 ● (!Cdc14 | (Cdc14:1 & Net1) | ((Cdc14:2 | 
Cdc14:3) & Net1:3)) & !Swi5 & !((Clb2 & !
 (Sic1 | Cdc6)) | (Clb5 & !Sic1) | Cln2 | ((Clb2 | 
Clb5) & (Cln3 | Bck2:2)) | (Clb5 & Clb2) | 
(Clb5:3 & Bck2))
● (!Cdc14 | (Cdc14:1 & Net1) | ((Cdc14:2 | 
Cdc14:3) & Net1:3)) & Swi5 & !((Clb2:3 & !
 (Sic1 | Cdc6)) | (((Clb2 & !(Sic1 | Cdc6)) | (Clb5 
& !Sic1)) & ((Cln2 & (Cln3 | Bck2)) | (Cln3 & 
Bck2))) | (((Clb2 & Clb5) | Clb2:3 | Clb5:3) & 
Cln2 & (Cln3 | Bck2)) | (Clb2:3 & Clb5:3))
● ((Cdc14:1 & !Net1) | (Cdc14:2 & (!Net1 | 
Net1:1))) & !Swi5 & !((Clb2 & !(Sic1 | Cdc6)) | 
(Clb5 & !Sic1) | Cln2)
● ((Cdc14:1 & !Net1) | (Cdc14:2 & (!Net1 | 
Net1:1))) & Swi5 & !(((Clb5 & Clb2 & !(Sic1 | 
Cdc6) & Cln2 & Cln3 & Bck2) | (Clb2:3 & !
 (Sic1 | Cdc6))) & (Clb5 | Cln2 | (Cln3 & Bck2)))
●Cdc14:3 & !Net1:3 & !Swi5 & !(Clb5 & Clb2 & 
Cln2 & Cln3 & Bck2)
●Cdc14:3 & !Net1:3 & Swi5 & !(Cln3 & Bck2 & 
Cln2 & Clb5 & !Sic1 & Clb2 & !(Sic1 | Cdc6))
Clb5 1 ●MASS & MBF_SBF & Cdc20
2 ●MASS & MBF_SBF & !Cdc20
Mpk1 1 ● !BUD
 Mih1 1 ● (Mpk1 & Clb2 & !(Sic1 | Cdc6)) | (!Mpk1 & (!
 Clb2 | Sic1 | Cdc6))
 2 ● !Mpk1 & Clb2 & !(Sic1 | Cdc6)
 Hsl1 1 ●BUD
 Swe1 1 ●MBF_SBF & ((Clb2 & !(Sic1 | Cdc6) & !Hsl1) | 
(Hsl1 & !Clb2) | Sic1 | Cdc6)
 2 ●MBF_SBF & !(Hsl1 | (Clb2 & !(Sic1 | Cdc6)))
 Clb2 1 ● ((MASS:1 & ((Swe1:1 & !Mih1:2) | (Swe1:2 & 
Mih1:1))) | (MASS & Swe1:2 & !Mih1)) & !
 Cdh1 & (!Cdc20 | (Cdc20:2 & Mcm1))
 2 ● ((MASS:1 & (!Swe1 | Mih1:2)) | (MASS:2 & (!
 Swe1:2 | Mih1))) & !Cdh1 & ((!Cdc20 & !
 Mcm1) | (Cdc20:2 & Mcm1))
3 ● ((MASS:1 & (!Swe1 | Mih1:2)) | (MASS:2 & (!
Swe1:2 | Mih1))) & !Cdh1 & !Cdc20 & Mcm1
Mcm1 1 ●Clb2 & !(Sic1 | Cdc6)
Mad2 1 ●ORI & !SPN
Cdc20 1 ● !Mad2 & Mcm1 & !(Clb2 & !(Sic1 | Cdc6))
 2 ● !Mad2 & Mcm1 & Clb2 & !(Sic1 | Cdc6)
 Cdc5Polo 1 ●Clb2 & !(Sic1 | Cdc6) & !Cdh1
 PP2ACdc55 1 ● !(!SeparaseEsp1 | (SeparaseEsp1:1 & 
SecurinPds1))
 2 ● !SeparaseEsp1 | (SeparaseEsp1:1 & 
SecurinPds1)
 Bub2Bfa1 1 ●ORI & (!SPN | PP2ACdc55 | !Cdc5Polo)
 Lte1 1 ● (Clb2 & !(Sic1 | Cdc6)) | SPN
 Tem1 1 ●Lte1
 Cdc15 1 ● !(Clb2 & !(Sic1 | Cdc6)) | (Cdc14 & !Net1)
 Net1 1 ● ((Cdc14 & Net1 & !(Clb2 & !(Sic1 | Cdc6)) & !
 PP2ACdc55) | (PP2ACdc55:1 & ((!Cdc14 & 
Clb2:2 & !(Sic1 | Cdc6)) | (Cdc14 & Clb2:2 & !
 (Sic1 | Cdc6))))) & !(((Cdc15:1 & Tem1) | 
Cdc15:2) & !Bub2Bfa1)
2 ● ((((Cdc14 & !Net1) | PP2ACdc55) & !(Clb2 & !
 (Sic1 | Cdc6))) | PP2ACdc55:2) & !(((Cdc15:1 & 
Tem1) | Cdc15:2) & !Bub2Bfa1)
Cdc14 1 ● (basal value)
Cdh1 1 ● (!Cdc14 | (Cdc14:1 & Net1) | ((Cdc14:2 | 
Cdc14:3) & Net1:3)) & !((Clb2 & !(Sic1 | Cdc6)) 
| (Clb5 & !Sic1) | (Cln3 & Cln2))
● ((Cdc14:1 & !Net1) | (Cdc14:2 & !Net1:3)) & !
 ((Clb5 & !Sic1 & Cln3 & ((Clb2 & !Cdc6) | 
Cln2)) | (Clb2:3 & !(Sic1 | Cdc6) & Cln3 & 
Cln2))
●Cdc14:3 & !Net1:3 & !((Clb5 & !Sic1 & Cln3 & 
Clb2 & !Cdc6 & Cln2) | (Clb2:3 & !(Sic1 | 
Cdc6) & Cln3 & Cln2))
BUD 1 ● (Cln2 | (Clb5 & !Sic1)) & !CYTOKINESIS
ORI 1 ● (Clb5 & !Sic1) | ((Clb2:2 | Clb2:3) & !(Sic1 | 
Cdc6)) | (ORI & (Clb5:3 | Clb2:3 | (Clb5:1 & !
 Sic1) | (Clb2:1 & !(Sic1 | Cdc6))))
SPN 1 ●Clb2 & !(Sic1 | Cdc6) & !CYTOKINESIS
SecurinPds1 1 ● (Mcm1 & MBF_SBF & !Cdc20) | ((Mcm1 | 
MBF_SBF) & !Cdh1 & !Cdc20)
SeparaseEsp1 1 ● (basal value)
MASS 2 ● !CYTOKINESIS
CYTOKINESIS 1 ●MASS & Clb2:2 & !(Sic1 | Cdc6)
2 ●MASS & ((Clb2:1 & CYTOKINESIS) | (Clb2:2 
& CYTOKINESIS & (Sic1 | Cdc6)))
""".split('\n')
while '' in text_coupled: text_coupled.remove('')

pmid = '19763337'

def standardize(text, file):
    full_text = []
    for line in text:
        if line[0] == ' ':
            line = line[1:]
        line = line.replace('●', ' ● ')
        if '●' in line:
            node = line.split(' ')
            while '' in node: node.remove('')
            if node[0] == '2' or node[0] == '3':
                state = node[0]
                node = full_text[len(full_text) - 1].split('●')[0].split(' ')
                while '' in node: node.remove('')
                node = node[0]
                full_text.append(node + ' ' + state + line[1:])
            elif node[0] == '●':
                node = full_text[len(full_text) - 1].split('●')[0]
                full_text.append(node + line)
    
            else:
                full_text.append(line)
        else:
            full_text[len(full_text) - 1] += line
    for i, line in enumerate(full_text):
        for j in range(1, 30):
            line = line.replace(''.join(' ' for k in range(30 - j)), ' ')
        full_text[i] = line
        if '(basal value)' in line:
            full_text[i] = ''
    while '' in full_text: full_text.remove('')
    
    r = []
    for i, line in enumerate(full_text):
        spline = line.split('●')
        node = spline[0].split(' ')
        while '' in node: node.remove('')
        node = node[0] + '=' + node[1]
        rule = spline[1].replace('(', ' ( ').replace(')', ' ) ').split(' ')
        while '' in rule: rule.remove('')
        for j, t in enumerate(rule):
            if t == '(' or t == ')':
                continue
            if t == '!':
                rule[j] = 'NOT'
                continue
            if t == '&':
                rule[j] = 'AND'
                continue
            if t == '|':
                rule[j] = 'OR'
                continue
            t = t.replace('!', 'NOT ')
            if ':' in t:
                t = t.replace(':', '=')
            else:
                t = 'NOT ' + t + '=0'
            rule[j] = t
        rule = ''.join(rule[j] + ' ' for j in range(len(rule))).replace('( ', '(').replace(' )', ')')
        r.append(node + ' :\t' + rule.replace('NOT NOT ', ''))
    
    g = open(pmid + file + '.txt', 'w')
    for line in r:
        g.write(line + '\n')
    g.close()

standardize(text_core, '_core_model')
standardize(text_morpho, '_morphogenesis_check_point')
standardize(text_exit, '_exit')
standardize(text_coupled, '')