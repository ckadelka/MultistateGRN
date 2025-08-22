#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 27 17:11:48 2025

@author: benco
"""

text = """
 TCR 1 APC
 CD28 1 APC
 IFNGR 1 IFNGR1 ∧ IFNGR2 ∧ (IFNG ∨ IFNGe)
 IFNGR1 1 (basalvalue)
 IFNGR2 1 (basalvalue)
 IL36R 1 IL36e ∧ IL1RL2 ∧ IL1RAP
 IL1RL2 1 (basalvalue)
 IL1R 1 (IL1Be ∨ IL1Ae) ∧ IL1RAP ∧ IL1R1
 IL1R1 1 STAT3
 IL1RAP 1 (basalvalue)
 IL2R 1 CGC ∧ IL2RB ∧ ¬IL2RA ∧ (IL2 ∨ IL2e)
 2 CGC ∧ IL2RB ∧ ¬IL2RA ∧ (IL2 ∨ IL2e)
 IL2RB 1 (basalvalue)
 IL2RA 1 (SMAD3 ∨ FOXP3 ∨ STAT5 ∨ NFKB) ∧ NFAT
 IL4R 1 CGC ∧ IL4RA ∧ (IL4 ∨ IL4e)
 IL6R 1 GP130 ∧ IL6RA ∧ (IL6e ∨ IL6)
 IL10R 1 IL10RA ∧ IL10RB ∧ (IL10 ∨ IL10e)
 IL12R 1 IL12RB1 ∧ IL12RB2 ∧ IL12e
 IL15R 1 CGC ∧ IL15RA ∧ IL2RB ∧ IL15e
 IL21R 1 GP130 ∧ CGC ∧ (IL21 ∨ IL21e)
 IL23R 1 GP130 ∧ IL12RB1 ∧ IL23e ∧ STAT3 ∧ RORGT
 IL27R 1 GP130 ∧ IL27RA ∧ IL27e
 IL27RA 1 (basalvalue)
 IFNAR 1 (IFNAe ∨ IFNBe) ∧ IFNAR1 ∧ IFNAR2
 IFNAR1 1 (basalvalue)
 IFNAR2 1 (basalvalue)
 TGFBR 1 TGFB ∨ TGFBe
 GP130 1 (basalvalue)
 IL6RA 1 (basalvalue)
 IL12RB1 1 (basalvalue)
 IL12RB2 1 ¬STAT6
 CGC 1 (basalvalue)
 IL10RA 1 (basalvalue)
 IL10RB 1 (basalvalue)
 IL4RA 1 (basalvalue)
 IL15RA 1 (basalvalue)
 IL29R 1 IL29e ∧ IL28RA ∧ IL10RB
 IL17RB 1 (basalvalue)
 IL18RAP 1 (basalvalue)
 IL18RA 1 (basalvalue)
 IL18R 1 IL18e ∧ IL18RAP ∧ IL18RA ∧ STAT4
 ST2 1 GATA3
 IL25R 1 IL17RB ∧ (IL25e ∨ IL25)
 IL33R 1 IL33e ∧ ST2 ∧ IL1RAP
 IL28RA 1 (basalvalue)
IKB 1 ¬TCR
NFKB 1 ¬IKB ∧ ¬FOXP3
NFAT 1 TCR ∧ CD28
TBET 1 (TBET ∨ STAT1 ∨ IL36R) ∧ ¬BCL6 ∧ ¬RORGT
GATA3 1 (¬GATA3 ∧ ¬TBET ∧ (STAT6 ∨ IL25R) ∧ ¬BCL6 ∧ ¬PU1 ∧ ¬IL29R) ∨ (GATA3 ∧ ¬BCL6 ∧ ¬PU1 ∧ ¬IL29R)
RORGT 1 TGFBR ∧ STAT3 ∧ ¬BCL6 ∧ ¬FOXP3
FOXP3 1 (STAT5 ∧ NFAT ∧ FOXP3 ∧ ¬STAT6) ∨ (STAT5 ∧ NFAT ∧ ¬FOXP3 ∧ SMAD3 ∧ ¬STAT1 ∧ ¬(STAT3 ∧ RORGT) ∧ ¬STAT6)
BCL6 1 ((STAT1 ∨ STAT3 ∨ STAT4) ∧ ¬TBET ∧ ¬STAT5) ∨ (STAT3 ∧ STAT4 ∧ ¬TBET)
STAT1 1 IFNAR ∨ IFNGR ∨ IL27R
STAT3 1 IL6R ∨ IL23R ∨ IL1R ∨ IL21R ∨ IL27R
STAT4 1 IL12R ∧ ¬GATA3
STAT5 1 ¬IL2R:2 ∧ (IL2R:1 ∨ IL15R)
2 IL2R:2
STAT6 1 IL4R
cMAF 1 TGFBR ∧ STAT3
PU1 1 TGFBR
SMAD3 1 TGFBR
IRF1 1 STAT1
RUNX3 1 TBET
IFNG 1 proliferation ∧ ¬FOXP3 ∧ NFAT ∧ ((TBET ∧ RUNX3) ∨ STAT4 ∨ IL18R)
IL4 1 NFAT ∧ proliferation ∧ GATA3 ∧ (STAT5 ∨ cMAF) ∧ ¬FOXP3 ∧ ¬((TBET ∧ RUNX3) ∨ IRF1)
IL2 1 (NFAT ∨ NFKB) ∧ ¬TBET ∧ ¬FOXP3 ∧ ¬(STAT5 ∧ STAT6)
IL17 1 NFAT ∧ proliferation ∧ RORGT ∧ NFKB ∧ STAT3 ∧ ¬(FOXP3 ∧ STAT1 ∧ STAT5 ∧ STAT6)
IL22 1 proliferation ∧ NFAT ∧ (STAT3 ∨ STAT1) ∧ ¬cMAF
IL9 1 (NFKB ∨ NFAT) ∧ proliferation ∧ (SMAD3 ∨ PU1 ∨ IL33R) ∧ STAT6
IL10 1 (GATA3 ∨ STAT3 ∨ STAT4 ∨ cMAF ∨ IRF1) ∧ NFAT ∧ proliferation ∧ ¬IL18R ∧ ¬IL33R
IL3 1 GATA3 ∧ proliferation ∧ NFAT
IL21 1 NFAT ∧ proliferation ∧ (STAT3 ∨ cMAF ∨ STAT4)
IL5 1 proliferation ∧ NFAT ∧ (GATA3 ∨ cMAF ∨ IL33R) ∧ ¬FOXP3
IL13 1 proliferation ∧ NFAT ∧ (GATA3 ∨ cMAF ∨ IL33R) ∧ ¬FOXP3
IL6 1 proliferation ∧ NFAT ∧ STAT3
TGFB 1 NFAT ∧ proliferation ∧ FOXP3
IL35 1 NFAT ∧ proliferation ∧ FOXP3
IL25 1 NFAT ∧ proliferation ∧ GATA3
IL31 1 NFAT ∧ proliferation ∧ STAT6
IL24 1 NFAT ∧ proliferation ∧ STAT6
proliferation 1 STAT5:2 ∨ proliferation
""".split('\n')
pmid= "25674559"

while '' in text: text.remove('')

for i, line in enumerate(text):
    if line == '':
        continue
    spline = line.split(' ')
    while '' in spline: spline.remove('')
    # the node identifier is not present, but must be same as the one for the previous line
    if spline[0] == '2':
        spline.insert(0, text[i - 1][0])
    text[i] = spline
r = []
binary_basal = []
for spline in text:
    rule = spline[0] + '=' + spline[1] + ' :\t'
    for i in range(2, len(spline)):
        t = spline[i].replace('¬', 'NOT ')
        if t == '(basalvalue)':
            rule = ''
            binary_basal.append(spline[0])
            break
        elif t == '∧':
            t = 'AND'
        elif t == '∨':
            t = 'OR'
        elif ':' in t:
            t = t.replace(':', '=')
        else:
            if t[len(t) - 1] == ')':
                t = t[:len(t) - 1] + '=1)'
            else:
                t += '=1'
        rule += t + ' '
    if not rule == '':
        r.append(rule)

for i in range(len(r)):
    for b in binary_basal:
        r[i] = r[i].replace(b + '=1', '')
    line = r[i].split(':\t')[1].split(' ')
    while '' in line: line.remove('')
    rule = ''
    for j in range(len(line)):
        term = line[j]
        if term == 'AND' or term == 'OR':
            if j == 0 or j == len(line) - 1 or line[j + 1] == 'AND' or line[j + 1] == 'OR' or 'NOT' in line[j - 1]:
                continue
            else:
                if rule == '':
                    continue
                rule += term + ' '
        else:
            rule += term + ' '
    r[i] = r[i].split(':\t')[0] + ':\t' + rule

g = open(pmid + '.txt', 'w')
for line in r:
    g.write(line + '\n')
g.close()