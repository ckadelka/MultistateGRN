#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 17:32:22 2025

@author: benco
"""

text = """
  EGFR1(EGFR_stimulus	
  |	
  SPRY)	
  &	
  !FGFR3	
  &	
  !GRB2(1)FGFR31!EGFR	
  &	
  FGFR3_stimulus	
  &	
  !GRB2(1)RAS1EGFR	
  |	
  FGFR3	
  |	
  GRB2(2)PTEN1TP53(3)PI3K1GRB2	
  &	
  RAS	
  &	
  !PTEN(3);	
  (4)AKT1PI3K(1)GRB21(FGFR3	
  &	
  !GRB2	
  &!SPRY)	
  |	
  EGFR(5)SPRY1RAS(5)E2F11!RB1	
  &	
  !RBL2	
  &((!(CHEK1_2:2&	
  ATM:2)	
  &	
  (RAS	
  |	
  E2F3))	
  |	
  (CHEK1_2:2&	
  ATM:2&	
  	
  !RAS	
  &E2F3:1))(6);	
  (7);	
  (8)2!RB1	
  &	
  !RBL2	
  &	
  ATM:2&	
  CHEK1_2:2&	
  (RAS	
  |	
  E2F3:2)(6);	
  (9)E2F31!RB1	
  &	
  !CHEK1_2:2&	
  RAS(10)2!RB1	
  &	
  CHEK1_2:2&	
  RAS(11)CyclinD11(RAS	
  |	
  AKT)	
  &	
  !p16INK4a	
  &	
  !p21CIP(12);	
  (13)CyclinE11!RBL2	
  &	
  !p21CIP	
  &	
  CDC25A	
  &	
  (E2F1	
  |	
  E2F3)(14);	
  (15)CyclinA1!RBL2	
  &	
  !p21CIP	
  &	
  CDC25A	
  &	
  (E2F1	
  |	
  E2F3)(14)CDC25A1!CHEK1_2	
  &	
  !RBL2	
  &	
  (E2F1	
  |	
  E2F3)(16);	
  (17)p16INK4a1Growth_inhibitors	
  &	
  !RB1(18);	
  (19)p14ARF1E2F1(20)RB11!CyclinD1	
  &	
  !CyclinE1	
  &	
  !p16INK4a	
  &	
  !CyclinA(21);	
  (22);	
  (18)RBL21!CyclinD1	
  &	
  !CyclinE1(23);	
  (24)p21CIP1!CyclinE1	
  &	
  (Growth_inhibitors|	
  TP53)	
  &	
  !AKT(15);	
  (25);	
  (26)ATM1DNA_damage	
  &	
  !E2F1(27)2DNA_damage	
  &	
  E2F1(28);	
  (29)CHEK1_21ATM	
  &	
  !E2F1(9);	
  (30)2ATM	
  &	
  E2F1(29)MDM21(TP53	
  |	
  AKT)	
  &	
  !p14ARF	
  &	
  !ATM	
  &	
  !RB1(31);	
  (32);	
  (33);	
  (34)TP531!MDM2	
  &	
  ((ATM	
  &CHEK1_2)	
  |	
  E2F1:2)(9);	
  (35);	
  (36);	
  (27)
  Proliferation1CyclinE1	
  |	
  CyclinA(26);	
  (37)Apoptosis1!E2F1:2&	
  TP53(38)2E2F1:2(39);	
  (27);	
  (40)Growth_arrest1p21CIP	
  |	
  RB1	
  |	
  RBL2
""".replace(';', '').replace('\n  ', ' ').replace('\t', ' ')
pmid = '26238783'

for i in range(1, 42):
    text = text.replace('(' + str(i) + ')', '#')
text = text.split('#')    

ntext = []
for line in text:
    if line[0] == '2':
        ntext.append(ntext[len(ntext) - 1].split('@')[0] + ' @2 ' + line[1:])
        continue
    for i in range(len(line)):
        if line[i] == '1':
            if line[i + 1] == '_' or line[i - 1] == 'p' or line[i - 2] == 'p' or line[i - 2] == '2' or line[i - 1] == 'D' or line[i - 1] == 'E' or line[i - 1] == 'B':
                continue
            else:
                ntext.append(line[0:i] + str(' @1 ') + line[i + 1:])
                break

for i in range(len(ntext)):
    newline = ntext[i].replace('&', ' AND ').replace('|', ' OR ').replace('!', ' NOT ').replace('\n', '').replace('\t', '')
    
    post = newline.split('@')[1][1:].split(' ')
    suffix = ''
    while '' in post: post.remove('')
    for x in post:
        if x.replace('(', '').replace(')', '') == '':
            suffix += ' ' + x + ' '
            continue
        if not ('OR' in x or 'AND' in x or 'NOT' in x):
            if ':' in x:
                x = x.replace(':', '=')
            elif x[len(x) - 1] == ')':
                if x[len(x) - 2] == ')':
                    x = x[:len(x) - 2] + '=1))'
                else:
                    x = x[:len(x) - 1] + '=1)'
            else:
                x = x + '=1'
        suffix += ' ' + x + ' '
    ntext[i] = (newline.split('@')[0].strip(' ') + '=' + newline.split('@')[1][0] + ' :\t' + suffix[1:]).replace('  ', ' ')

g = open(pmid + '.txt', 'w')
for line in ntext:
    g.write(line + '\n')
g.close()