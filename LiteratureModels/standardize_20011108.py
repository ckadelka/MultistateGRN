# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 14:00:11 2025

@author: benco
"""

text = """
housekeeping = FADD 0 Harper et al, 2003, Muzio 
et al, 1998 Homo sapiens 
housekeeping = IKKdeact 0 Devin et al, 2000 Homo sapiens,  
Mus musculus 
housekeeping = IRS 0 reviewed by White 1998 Homo sapiens, 
Mus musculus 
housekeeping = P14-3-3 0 Masters et al, 2001 Homo sapiens 
housekeeping = proC8 0 Lavrik et al, 2003 Homo sapiens 
housekeeping = RIP 0 Wang et al, 2008, Devin et 
al, 2000 
Homo sapiens,  
Mus musculus 
housekeeping = TRADD 0 Hsu et al, 1995, Hsu et al, 
1996 Homo sapiens 
comp1 + IKKdeact = comp1-IKK* 2 Devin et al, 2000 Homo sapiens,  
Mus musculus 
comp1 = NIK 2 Senftleben et al, 2001, 
Xiao et al, 2001 
Homo sapiens, 
Mus musculus 
TNF = JNK 2 Liu et al, 1996 Homo sapiens 
TNF = TNFR-1 2 Tartaglia et al, 1993 Homo sapiens 
TNFR-1 + RIP + TRADD + TRAF-2 + c-IAP = comp1 2 
Boldin et al, 1996, Hsu et 
al, 1995, Rothe et al, 1995, 
Mahoney et al, 2008 
Homo sapiens, 
Mus musculus 
C8*-comp2= C8* 3 Wang et al, 2008, Chen 
and Goeddel, 2002 
Homo sapiens, 
Mus musculus 
proC8 + comp2 = C8*-comp2 3 Chen and Goeddel, 2002 Homo sapiens, 
Mus musculus 
RIP-deubi + comp1 + FADD = comp2 3 Boldin et al, 1996, Wang et 
al, 2008 
Homo sapiens, 
Mus musculus 
smac = RIP-deubi 3 Wang et al, 2008 Mus musculus 
smac-mimetics = smac 3 Li et al, 2004, Wang et al, 
2008 Homo sapiens 
!BAD-14-3-3 + !Bcl-xl = Bax 4 Cheng et al, 2001 Mus musculus 
!I-kBa + !I-kBe = NF-κB 4 Baeuerle and Baltimore, 
1988, Kearns et al, 2006 
Homo sapiens, 
Mus musculus 
!PARP + !ICAD = CAD 4 
Enari et al, 1998, Sakahira 
et al, 1998, reviewed by 
Cryns and Yuan, 1998 
Homo sapiens, 
Mus musculus 
!PKB = GSK-3 4 
Datta et al, 1997, reviewed 
by Nystrom and Quon, 
1999 
Homo sapiens, 
Mus musculus 
!UV = Bcl-xl 4 Zhang and Rosdahl, 2005 Homo sapiens 
1 C3*p20 = 1 C3*p17 4 Han et al, 1997 Homo sapiens 
1 C8* = 1 C3*p20 4 Han et al, 1997 Homo sapiens 
1 C8*-DISC* = 1 C8* 4 Lavrik et al, 2003 Homo sapiens 
1 DISC* + proC8 = 1 C8*-DISC* 4 Lavrik et al, 2003 Homo sapiens 
1 Fas + FADD =1 DISC* 4 Boldin et al, 1996 Homo sapiens 
1 FasL = 1 Fas, 2 FasL = 2 Fas 4 Lavrik et al, 2007, Itoh and 
Nagata, 1993 
Homo sapiens, 
Mus musculus 
2 C3*p20 + 2 !XIAP = 2 C3*p17 4 
Deveraux et al, 1997, Li et 
al, 2004, Chai et al, 2000, 
Verhagen et al, 2000, 
Rehm et al, 2006 
Homo sapiens, 
Mus musculus 
2 C8* = 2 C3*p20 4 Han et al, 1997 Homo sapiens 
 - 5 - 
2 C8*-DISC* = 2 C8* 4 Lavrik et al, 2003 Homo sapiens 
2 DISC* + proC8 = 2 C8*-DISC* 4 Lavrik et al, 2003 Homo sapiens 
2 FAS + FADD = 2 DISC* 4 Boldin et al, 1996 Homo sapiens 
AdCy = cAMP 4 Jelinek et al, 1993 
Homo sapiens, 
Rattus 
norvegicus 
Apaf-1 = apopto 4 Zhou et al, 1997, Deveraux 
et al, 1998, Li et al, 1997 Homo sapiens 
apopto = C9* 4 
Zhou et al, 1997, Li et al, 
1997, Slee et al, 1999, 
Srinivasula et al, 1998 
Homo sapiens 
BAD + P14-3-3 = BAD-14-3-3 4 Datta et al, 1997 Homo sapiens, 
Mus musculus 
Bax = cyt-c 4 
Madesh et al, 2002, Wang 
et al, 1996, Desagher et al, 
1999 
Homo sapiens 
Bax = smac 4 
Madesh et al, 2002, Wang 
et al, 1996, Desagher et al, 
1999 
Homo sapiens, 
Mus musculus 
C3*p17 + c-IAP = C3*-c-IAP 4 Deveraux et al, 1997 Homo sapiens 
C3*p17 = C3*-XIAP 4 Deveraux et al, 1997 Homo sapiens 
C3*-XIAP = BIR1-2 4 Rehm et al, 2006 Homo sapiens 
C8* + !XIAP = C3*p17 4 
Li et al, 2004, Chai et al, 
2000, Verhagen et al, 
2000, Rehm et al, 2006 
Homo sapiens, 
Mus musculus 
C8* + P = tBid 4 Zhao et al, 2003, Roucou 
et al, 2002 Mus musculus 
C8*-DISC* + !FLIP = C8* 4 
Lavrik et al, 2007, Irmler 
et al, 1997, Tschopp et al, 
1998 
Homo sapiens, 
Mus musculus 
C8*-DISC* + FLIP = C8*-FLIP 4 
Lavrik et al, 2007, Irmler 
et al, 1997, Tschopp et al, 
1998 
Homo sapiens 
C9* = 2 C3*p20 4 Slee et al, 1999 Homo sapiens 
CAD + !PARP + !gelsolin = apoptosis 4 
Enari et al, 1998, Sakahira 
et al, 1998, reviewed by 
Cryns and Yuan 1998 
Homo sapiens, 
Mus musculus 
cAMP = PKA 4 Cheng et al, 1998 Homo sapiens, 
Mus musculus 
cyt-c + smac = Apaf-1 4 Zhou et al, 1997 Homo sapiens 
Glucagon = GR 4 Jelinek et al, 1993 
Homo sapiens, 
Rattus 
norvegicus 
GR = AdCy 4 Jelinek et al, 1993 
Homo sapiens, 
Rattus 
norvegicus 
Grb2-SOS = Ras 4 reviewed by Kolch 2000 Homo sapiens, 
Mus musculus 
housekeeping + !C3*p17 = gelsolin 4 Kothakota et al, 1997, 
Janicke et al, 1998 
Homo sapiens, 
Mus musculus 
housekeeping + !IKK* + !comp1-IKK* +  
2 !UV = I-kBa 4 Mercurio et al, 1997, 
DiDonato et al, 1997 Homo sapiens 
housekeeping + !IKK* + !comp1-IKK* +  
2 !UV = I-kBb 4 Mercurio et al, 1997, 
DiDonato et al, 1997 Homo sapiens 
 - 6 - 
housekeeping + !IKK* + !comp1-IKK* +  
2 !UV = I-kBe 4 Mercurio et al, 1997, 
DiDonato et al, 1997 
Homo sapiens 
 
housekeeping + !P = Bid 4 McKenzie et al, 2008, 
Krams et al, 2007 Mus musculus 
housekeeping + 2 !C3*p17 = PARP 4 
Enari et al, 1998, Lazebnik 
et al, 1994, Sakahira et al, 
1998 
Homo sapiens, 
Mus musculus 
housekeeping + 2!C3*p17 = ICAD 4 Kothakota et al, 1997, 
Janicke et al, 1998 
Homo sapiens, 
Mus musculus 
IL-1 + IKKdeact = IKK* 4 
Poeppelmann et al, 2005, 
Barisic et al, 2008, 
Kothny-Wilkes et al, 1999 
Homo sapiens 
Insulin + !IRS-P2 = IR 4 reviewed by Gual et al, 
2005 
Homo sapiens, 
Mus musculus 
IR + IRS = IRS-P 4 
reviewed by White 1998, 
Saltiel and Kahn, 2001, 
reviewed by Nystrom and 
Quon, 1999, reviewed by 
White 2002 
Homo sapiens, 
Mus musculus 
IRS-P = PI3K 4 
Saltiel and Kahn, 2001, 
reviewed by Nystrom and 
Quon, 1999, reviewed by 
White 2002 
Homo sapiens, 
Mus musculus 
IRS-P = Shc 4 Saltiel and Kahn, 2001, 
reviewed by White 2002 
Homo sapiens, 
Mus musculus 
NIK + IKKdeact = IKK* 4 Senftleben et al, 2001, 
Xiao et al, 2001 
Homo sapiens, 
Mus musculus 
Pak1 = MEK 4 reviewed by Kolch 2000 Homo sapiens, 
Mus musculus 
Pak1 = Raf 4 reviewed by Kolch 2000 Homo sapiens, 
Mus musculus 
PDK1 = PKB 4 Good et al, 1998, 
Vanhaesebrock et al, 2000 Homo sapiens 
PDK1 = PKC 4 Pruett et al, 1995, Good et 
al, 1998 Homo sapiens 
PI3K + IKKdeact = IKK* 4 Datta et al, 1997 Homo sapiens, 
Mus musculus 
PI3K = PIP3 4 
Saltiel and Kahn, 2001, 
reviewed by Nystrom and 
Quon, 1999, reviewed by 
White 2002 
Homo sapiens, 
Mus musculus 
PI3K = Rac 4 
Saltiel and Kahn, 2001, 
reviewed by Nystrom and 
Quon, 1999, reviewed by 
White 2002 
Homo sapiens, 
Mus musculus 
PIP3 = PDK1 4 
Vanhaesebrock et al, 2000, 
Saltiel and Kahn, 2001, 
reviewed by Nystrom and 
Quon, 1999, reviewed by 
White 2002 
Homo sapiens, 
Mus musculus 
PKB = BAD 4 Datta et al, 1997 Homo sapiens, 
Mus musculus 
PKC = Raf 4 reviewed by Kolch 2000 Homo sapiens, 
Mus musculus 
Rac = Pak1 4 
Saltiel and Kahn, 2001, 
reviewed by Nystrom and 
Quon, 1999, reviewed by 
White 2002 
Homo sapiens, 
Mus musculus 
Raf + IKKdeact = IKK* 4 Lin et al, 1999 Homo sapiens 
 - 7 - 
Ras = ERK1/2 4 reviewed by Kolch 2000 Homo sapiens, 
Mus musculus 
Ras = p38 4 reviewed by Kolch 2000 Homo sapiens, 
Mus musculus 
Shc = Grb2-SOS 4 
Pruett et al, 1995, Saltiel 
and Kahn, 2001, reviewed 
by White 2002 
Homo sapiens, 
Mus musculus 
smac = smac-XIAP 4 
Li et al, 2004, Chai et al, 
2000, Verhagen et al, 
2000, Rehm et al, 2006 
Homo sapiens, 
Mus musculus 
smac-XIAP + 1 C8* = 2 C3*p17 4 
Li et al, 2004, Chai et al, 
2000, Verhagen et al, 
2000, Rehm et al, 2006 
Homo sapiens, 
Mus musculus 
T2R + C8 = C8* 4 based on Scaffidi et al, 
1998, Walter et al, 2008 
Homo sapiens, 
Mus musculus 
T2R = P 4 based on Scaffidi et al, 
1998, Walter et al, 2008 
Homo sapiens, 
Mus musculus 
T2RL = T2R 4 based on Scaffidi et al, 
1998, Walter et al, 2008 
Homo sapiens, 
Mus musculus 
tBid = Bax 4 
Wei et al, 2001, Madesh et 
al, 2002, Wang et al, 1996, 
Desagher et al, 1999, Zhao 
et al, 2003 
Homo sapiens, 
Mus musculus 
TNF + IRS = IRS-P2 4 Gual et al, 2005 Homo sapiens, 
Mus musculus 
UV = Bax 4 Chen et al, 2007 Mus musculus 
!BIR1-2 = C3*p17 5 Rehm et al, 2006, 
Deveraux et al, 1997 Homo sapiens 
!C3*p17 = Bcl-xl 5 Cheng et al, 1998, Clem et 
al, 1998 
Homo sapiens, 
Rattus 
norvegicus 
!PKA = Raf 5 reviewed by Kolch 2000 Homo sapiens, 
Mus musculus 
!PKC = IRS-P 5 Fea and Roth, 1997, 
Ravichandran et al, 2001 
Homo sapiens, 
Mus musculus 
!smac = c-IAP 5 Wang et al, 2008 Homo sapiens 
2 C3*p17 = C9* 5 
Zhou et al, 2003, Rehm et 
al, 2006, Deveraux et al, 
1997 
Homo sapiens 
C3*p17 = C6 5 Cowling and Downward 
2002, Murphy et al, 2004 Homo sapiens 
C6 = C8* 5 Murphy et al, 2004 Homo sapiens 
PKB = IRS-P 5 
Paz et al, 1999, reviewed 
by Nystrom and Quon, 
1999 
Homo sapiens, 
Rattus 
norvegicus 
!A20 = comp1-IKK* 10 Song et al, 1996, Wertz et 
al, 2004, Lee et al,2000 
Homo sapiens, 
Mus musculus 
!NF-κB = JNK 10 Tang et al, 2001 Mus musculus 
housekeeping = c-IAP 10 Poeppelmann et al, 2005, 
Barisic et al, 2008 Homo sapiens 
housekeeping = FLIP 10 Poeppelmann et al, 2005 Homo sapiens 
housekeeping = XIAP 10 Poeppelmann et al, 2005 Homo sapiens 
housekeeping = TRAF-2 10 Krikos et al, 1992 Homo sapiens 
NF-κB = 2 XIAP 10 Barisic et al, 2008 Homo sapiens 
NF-κB = 2 c-IAP 10 Barisic et al, 2008 Homo sapiens 
 - 8 - 
NF-κB = 2 FLIP 10 
Kreuz et al, 2001, Micheau 
et al, 2001, Barisic et al, 
2008 
Homo sapiens, 
Mus musculus 
NF-κB = A20 10 Song et al, 1996, Krikos et 
al, 1992 
Homo sapiens,  
Mus musculus 
NF-κB = Bcl-xl 10 
Lee et al, 1999, Tamatani 
et al, 1999, Chen et al, 
2000 
Homo sapiens, 
Rattus 
norvegicus 
NF-κB = I-kBa  10 Tian et al, 2005 Homo sapiens 
NF-κB = I-kBb 10 Tian et al, 2005 Homo sapiens 
NF-κB = I-kBe  10 Rothe et al, 1995, Wang et 
al, 1998 
Homo sapiens, 
Mus musculus 
NF-κB = TRAF-2 10 
Poeppelmann et al, 2005, 
Barisic et al, 2008, 
Kothny-Wilkes et al, 1999 
Homo sapiens
""".split('\n')
while '' in text: text.remove('')

text_add_1 = """
!A20 = comp1-IKK* 
!BIR1-2 = C3*p17 
!NF-κB = JNK 
!PKA = Raf 
!PKC = IRS-P 
!smac = c-IAP 
2 C3*p17 = C9* 
C3*p17 = C6 
C6 = C8* 
NF-κB = I-kBa 
NF-κB = I-kBb 
NF-κB = I-kBe 
PKB = IRS-P
""".split('\n')
while '' in text_add_1: text_add_1.remove('')

text_add_2 = """
1 C3*p20 = 1 C3*p17 
1 C8* = 1 C3*p20 
1 C8*-DISC* = 1 C8* 
1 DISC* + proC8 = 1 C8*-DISC* 
1 FAS + FADD = 1 DISC* 
1 FASL = 1 FAS 
2 C3*p17 = C9* 
2 C3*p20 + 2 !XIAP = 2 C3*p17 
2 C8* = 2 C3*p20 
2 C8*-DISC* = 2 C8* 
2 DISC* + proC8 = 2 C8*-DISC* 
2 FAS + FADD = 2 DISC* 
2 FASL = 2 FAS 
smac-XIAP + 1 C8* = 2 C3*p17
""".split('\n')
while '' in text_add_2: text_add_2.remove('')

pmid = '20011108'

# note: Regulators are on the left
culled_text = []
for line in text:
    if '=' in line:
        eq = line.split('=')
        right = eq[1].split(' ')
        if right[2].isdigit():
            right = ''.join(' ' + right[i] for i in range(0, 2))
        elif right[3].isdigit():
            right = ''.join(' ' + right[i] for i in range(0, 3))
        culled_text.append(' ' + (eq[0] + ' = ' + right).replace('   ', ' ').replace('  ', ' '))
for line in text_add_1:
    culled_text.append(' ' + line)
for line in text_add_2:
    culled_text.append(' ' + line)
for i, line in enumerate(culled_text):
    culled_text[i] = line.replace(' 1 ', ' 1:').replace(' 2 ', ' 2:').replace(' 2!', ' 2:!')

rdict = dict()
for i, line in enumerate(culled_text):
    node = line.split('=')[1].replace(' ', '')
    if ':' in node:
        node = node.split(':')[1] + '=' + node.split(':')[0]
    else:
        node = node + '=1'    
    rule = line.split('=')[0][1:].replace('+', 'AND')
    rs = rule.split(' ')
    while '' in rs: rs.remove('')
    for j, r in enumerate(rs):
        if r == 'AND':
            continue
        if '!' in r:
            if ':' in r:
                r = r.replace('!', '')
                r = 'NOT ' + r.split(':')[1] + ':' + r.split(':')[0]
            else:
                r = r.replace('!', '') + ':0'
        else:
            if ':' in r:
                r = r.split(':')[1] + ':' + r.split(':')[0]
            else:
                r = r + ':1'
        rs[j] = r
    rule = ''.join(rs[j] + ' ' for j in range(len(rs))).replace(':', '=')
    
    try:
        rdict[node] = '(' + rdict[node] + ') OR (' + rule + ')'
    except KeyError:
        rdict[node] = rule

sort = sorted(rdict)
g = open(pmid + '.txt', 'w')
for key in sort:
    g.write(key + ' :\t' + rdict[key] + '\n')
g.close()