#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 27 18:29:26 2025

@author: benco
"""

import numpy as np

text = """
x1= LIP
 x2+2x2²x23+2x2x23²+2x2²x23²+2x2²x3+x2²x23x3+2x2²x23²x3+2x2x3²+
 2x2²x3²+2x2²x23x3²+x2x23²x3²+2x2²x4+x2²x23x4+2x2²x23²x4+x2²x3x4+
 2x2²x23x3x4+x2²x23²x3x4+2x2²x3²x4+x2²x23x3²x4+2x2²x23²x3²x4+2x2x4²+
 2x2²x4²+2x2²x23x4²+x2x23²x4²+2x2²x3x4²+x2²x23x3x4²+2x2²x23²x3x4²+
 x2x3²x4²+2x2²x23x3²x4²+2x2x23²x3²x4²+2x2²x23²x3²x4²+x8+2x2x8+x2²x8+
 2x23²x8+x2x23²x8+2x2²x23²x8+2x3²x8+x2x3²x8+2x2²x3²x8+x23²x3²x8+
 2x2x23²x3²x8+x2²x23²x3²x8+2x4²x8+x2x4²x8+2x2²x4²x8+x23²x4²x8+
 2x2x23²x4²x8+x2²x23²x4²x8+x3²x4²x8+2x2x3²x4²x8+x2²x3²x4²x8+
 2x23²x3²x4²x8+x2x23²x3²x4²x8+2x2²x23²x3²x4²x8+x2x8²+x2²x8²+2x23x8²+
 x2²x23x8²+2x23²x8²+2x2x23²x8²+2x3x8²+x2²x3x8²+x23x3x8²+2x2²x23x3x8²+
 2x23²x3x8²+x2²x23²x3x8²+2x3²x8²+2x2x3²x8²+2x23x3²x8²+x2²x23x3²x8²+
 x2x23²x3²x8²+x2²x23²x3²x8²+2x4x8²+x2²x4x8²+x23x4x8²+2x2²x23x4x8²+
 2x23²x4x8²+x2²x23²x4x8²+x3x4x8²+2x2²x3x4x8²+2x23x3x4x8²+x2²x23x3x4x8²+
 x23²x3x4x8²+2x2²x23²x3x4x8²+2x3²x4x8²+x2²x3²x4x8²+x23x3²x4x8²+
 2x2²x23x3²x4x8²+2x23²x3²x4x8²+x2²x23²x3²x4x8²+2x4²x8²+2x2x4²x8²+
 2x23x4²x8²+x2²x23x4²x8²+x2x23²x4²x8²+x2²x23²x4²x8²+2x3x4²x8²+x2²x3x4²x8²+
 x23x3x4²x8²+2x2²x23x3x4²x8²+2x23²x3x4²x8²+x2²x23²x3x4²x8²+x2x3²x4²x8²+
 x2²x3²x4²x8²+2x23x3²x4²x8²+x2²x23x3²x4²x8²+2x23²x3²x4²x8²+2x2x23²x3²x4²x8²
x2= TfR1
 x19²+x2+2x19²x2+2x2²+x19x2²+x5²+2x19²x5²+2x2x5²+
 x19²x2x5²+x2²x5²+2x19²x2²x5²+x2²x6+2x19x2²x6+x19²x2²x6+
 x6²+2x19²x6²+2x2x6²+x19²x2x6²+x19x2²x6²+x19²x2²x6²+2x5²x6²+
 x19²x5²x6²+x2x5²x6²+2x19²x2x5²x6²+2x2²x5²x6²+x19²x2²x5²x6²
x3= Fpn
 1+x3²+2x3²x5²+2x6+x3x6+2x3²x6+x6²+2x3x6²+x3²x5²x6²+2x7+x3x7+2x3²x7+
 x6x7+2x3x6x7+x3²x6x7+2x6²x7+x3x6²x7+2x3²x6²x7+x7²+2x3x7²+x3²x5²x7²+
 2x6x7²+x3x6x7²+2x3²x6x7²+x6²x7²+2x3x6²x7²+2x3²x6²x7²+2x3²x5²x6²x7²
x4= Ft
 1+x4²+2x4²x5²+2x6+x4x6+2x4²x6+x6²+2x4x6²+x4²x5²x6²
x5= IRP1
 1+2x1+x1²+x1x5+2x1²x5+x5²+2x1x5²
x6= IRP2
 1+2x1+x1²+x1x19²+2x1²x19²+x1x6+2x1²x6+2x1x19²x6+
 x1²x19²x6+x6²+2x1x6²+x1²x19x6²+x1x19²x6²+x1²x19²x6²
x7= Hep
 x15²+x7+2x15²x7+2x7²+x15x7²
x8= HO-1
 x10²+x13²+2x10²x13²+x8+2x10²x8+2x13²x8+x10²x13²x8+2x8²+
 x10x8²+x13x8²+2x10x13x8²+x10²x13x8²+x10x13²x8²+x10²x13²+x8²
x9= ALAS1
 x22²+2x10x22²+x10²x22²+x9+2x22²x9+x10x22²x9+
 2x10²x22²x9+2x9²+x22x9²+2x10²x22x9²+2x10x22²x9²+2x10²x22²x9²
x10= Heme
 x9+2x8²x9+2x8x9²+2x8²x9²
x11= ROS
 x1+2x1²x12+2x1x12²+2x1²x12²+x16+2x1x16+x1²x16+2x12²x16+x1x12²x16+
 2x1²x12²x16+x1x16²+x1²x16²+2x12x16²+x1²x12x16²+2x12²x16²+2x1x12²x16²+
 x21+2x1x21+x1²x21+2x12²x21+x1x12²x21+2x1²x12²x21+2x16x21+x1x16x21+
 2x1²x16x21+x12²x16x21+2x1x12²x16x21+x1²x12²x16x21+x16²x21+2x1x16²x21+
 x1²x16²x21+2x12²x16²x21+x1x12²x16²x21+2x1²x12²x16²x21+x1x21²+x1²x21²+
 2x12x21²+x1²x12x21²+2x12²x21²+2x1x12²x21²+x16x21²+2x1x16x21²+
 x1²x16x21²+2x12²x16x21²+x1x12²x16x21²+2x1²x12²x16x21²+x16²x21²+
 x1x16²x21²+x12x16²x21²+2x1²x12x16²x21²+2x1x12²x16²x21²+2x1²x12²x16²x21²
x12= Antioxidant Enzymes (AE)
 x12+2x12²+x12²x13+x13²+2x12x13²
x13= Nrf2
 1+x13²+2x14+x13x14+2x13²x14+x14²+2x13x14²+x13²x14²x16+x14x16²+
 2x13x14x16²+x13²x14x16²+2x14²x16²+x13x14²x16²+x13²x14²x16²+x13²x14²x18+
 2x13²x14²x16x18+x13²x14²x16²x18+x14x18²+2x13x14x18²+x13²x14x18²+
 2x14²x18²+x13x14²x18²+x13²x14²x18²+x13²x14²x16x18²+2x14x16²x18²+
 x13x14x16²x18²+2x13²x14x16²x18²+x14²x16²x18²+2x13x14²x16²x18²
x14= Keap1
 1+2x11+x11²+x11x14+2x11²x14+2x11x14²+
 x11²x14²+x13x14²+2x11²x13x14²+2x13²x14²+x11²x13²x14²
x15= IL-6
 1+x15²+2x8+x11²x8+x15x8+2x11²x15x8+2x15²x8+x11²x15²x8+
 x8²+2x11²x8²+2x15x8²+x11²x15x8²+x11x15²x8²+x11²x15²x8²
x16= Ras
 x15²+x16+2x15²x16+2x16²+x15x16²+x16²x17+2x15x16²x17+x15²x16²x17+
 x17²+2x15²x17²+2x16x17²+x15²x16x17²+x15x16²x17²+x15²x16²x17²+
 2x15²x20+x15²x16x20+2x15²x16²x20+2x17²x20+x15²x17²x20+
 x16x17²x20+2x15²x16x17²x20+2x16²x17²x20+x15²x16²x17²x20+
 x15²x20²+2x15²x16x20²+2x15x16²x20²+2x15²x16²x20²+2x16²x17x20+
 x15x16²x17x20²+2x15²x16²x17x20²+x17²x20²+2x15²x17²x20²+
 2x16x17²x20²+x15²x16x17²x20²+2x16²x17²x20²+2x15x16²x17²x20²
x17= SOS
 1+x17²+2x18+x17x18+2x17²x18+x18²+2x17x18²+x17²x18²x21+
 x18x21²+2x17x18x21²+x17²x18x21²+2x18²x21²+x17x18²x21²+x17²x18²x21²
x18= ERK
 x16²+x18+2x16²x18+2x18²+x16x18²
x19= c-Myc
 x18²+x19+2x18²x19+2x19²+x18x19²
x20= GAPs
 x20+2x20²+x20²x21+x21²+2x20x21²
x21= EGFR
 x11²+x21+2x11²x21+2x21²+x11x21²
x22= LIPmt
 x23+2x10²x23+2x10x23²+2x10²x23²+2x23²x24+x10x23²x24+
 2x10²x23²x24+2x23x24²+x10²x23x24²+2x23²x24²+2x10x23²x24²
x23= Mfrn
 1+2x22+x22²+x22x23+2x22²x23+x23²+2x22x23²
x24= Ftmt
 x22²+x24+2x22²x24+2x24²+x22x24²
""".replace('Antioxidant Enzymes (AE)', 'AE').split('\n')
pmid = "28166223"

while '' in text: text.remove('')
header_lines = []
for i, line in enumerate(text):
    if line == '':
        continue
    if '=' in line:
        header_lines.append(i)
header_lines.append(len(text))

name_dict = dict()
r = []
for i in range(len(header_lines) - 1):
    header_split = text[header_lines[i]].split('=')
    name_dict[header_split[0].replace('x', '')] = header_split[1][1:] # remove leading space from name
    
    j  = header_lines[i] + 1
    j_esc = header_lines[i + 1]
    polynomial = ''
    while j < j_esc:
        polynomial += text[j].replace(' ', '')
        j += 1
    poly_split = polynomial.split('+')
    
    for k, term in enumerate(poly_split):
        if 'x' in term:
            term_split = term.split('x')
            while '' in term_split: term_split.remove('')
            reformed_term = ''
            for l, ts in enumerate(term_split):
                if l == 0 and not term[0] == 'x':
                    reformed_term += ts
                    continue
                if l > 0:
                    reformed_term += ' * '
                reformed_term += 'x[' + ts + ']'
            reformed_term = reformed_term.split(' ')
            for l, ts in enumerate(reformed_term):
                if '²' in ts:
                    ts = ts.replace('²', '')
                    ts += ' * ' + ts
                    reformed_term[l] = ts
            poly_split[k] = ''.join(reformed_term[l] + ' ' for l in range(len(reformed_term)))
    
    r.append(''.join(poly_split[k] + '+ ' for k in range(len(poly_split)))[:-3]) # remove trailing ' + '

# Note that every node in this model is represented as a ternary
rule_eval = []
for i, rule in enumerate(r):
    num_set = []
    rsplit = rule.split('[')
    for j, sr in enumerate(rsplit):
        if ']' in sr:
            num_set.append(sr.split(']')[0])
    num_set = set(num_set)
    for j in range(3**len(num_set)):
        trinary = ''
        decimal = j
        while decimal > 0:
            trinary += str(decimal % 3)
            decimal //= 3
        trinary = ("{:0" + str(len(num_set)) + "s}").format(trinary) # pad with zeros
        x = np.zeros(25, int)
        formatted_rule = ''
        for k, num in enumerate(num_set):
            x[int(num)] = int(trinary[k])
            formatted_rule += name_dict[num] + '=' + trinary[k] + ' AND '
        rule_eval.append(name_dict[str(i + 1)] + '=' + str(eval(rule) % 3) + ' :\t' + formatted_rule[:-5]) # remove trailing ' AND '

# Convert to DNF, since instead of deconstructing the polynomial to determine the logical
# equivalence, we are simply evaluating it for every possible value, and thereby get a
# truth table-esque result, and the simplest way to combine the assortment of conjunctions
# is through their conjunction
# Do note, however, that this is not optimal space-wise nor time-wise, and it may be preferable
# to either rework this in the future to actually derive the logical equivalence in the
# polynomial or to simplify the DNF representation via an analogous 'combining of like terms'
rdict = dict()
for rule in rule_eval:
    rsplit = rule.split(':\t')
    try:
        if 'OR' in rdict[rsplit[0]]:
            rdict[rsplit[0]] = rdict[rsplit[0]] + ' OR (' + rsplit[1] + ')'
        else:
            rdict[rsplit[0]] = '(' + rdict[rsplit[0]] + ') OR (' + rsplit[1] + ')'
    except KeyError:
        rdict[rsplit[0]] = rsplit[1]

# Sort, so that each level of activation (0, 1, 2) for a node is contiguously
# listed instead of being sporatically scattered all over the place
sort_order = sorted(rdict)

g = open(pmid + '.txt', 'w')
for key in sort_order:
    g.write(key + ':\t' + rdict[key] + '\n')
g.close()