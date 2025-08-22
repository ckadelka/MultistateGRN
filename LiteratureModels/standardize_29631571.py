# -*- coding: utf-8 -*-
"""
Created on Thu Apr  3 15:28:47 2025

@author: benco
"""

text = """
Injury (2)* = 0 if M2 = 1 and previous injury = 1
DAMPs (3) = 0 if Injury = 0 AND ROS = 0 regardless of M2
DAMPs =0 if (Injury = 1 XOR** ROS = 1) and M2 = 1
DAMPs =1 if (Injury = 1 XOR** ROS = 1) and M2 = 0 unless previous DAMPs = 2
DAMPs =1 if (Injury = 1 AND ROS = 1) and M2 = 1
DAMPs =2 if (Injury = 1 AND ROS = 1) and M2 = 0
M1 (3) = 0 if (CCL2 = 0)
M1 = 1 if CCL2 = 1
M1 = 2 CCL2 = 2
M2 (2) = 1 if M1 = 1
M2 = 0 otherwise
CD13 (2)* = 1 if DAMPs = 1 or 2
CD13 = 0 otherwise
TRIF (3) = 0 if DAMPs = 0 regardless of CD13
TRIF = 1 if (DAMPs = 1) and (CD13 = 1)
TRIF = 2 if (DAMPs = 1) and (CD13 = 0)
TRIF = 2 if DAMPs = 2 regardless of CD13
IRF3 (3) = TRIF (3)
IFN-β (3) = IRF3
ROS (2) = 1 IFNβ = 2
ROS = 0 otherwise
MyD88 = DAMPs (3)
pIRAK = MyD88 (3)
NF-kB = 0 if M2 = 1 and (pIRAK = 0 or 1)
NF-kB = pIRAK (3) otherwise
CCL2 = NF-kB (3)
""".replace('previous', '').replace('\u2009', ' ').replace('*', '').split('\n')
pmid = '29631571'

while '' in text: text.remove('')

rdict = dict()
for line in text:
    pre = line.split('=')[0].replace(' ', '')
    if '(2)' in pre or '(3)' in pre:
        states = int(pre.split('(')[1].split(')')[0])
        for i in range(states):
            rdict[pre.split('(')[0] + '=' + str(i)] = ''
    post = line.split('=')[1].replace(' ', '')
    if '(2)' in post or '(3)' in post:
        states = int(post.split('(')[1].split(')')[0])
        for i in range(states):
            rdict[pre.split('(')[0] + '=' + str(i)] = ''

for line in text:
    node = line.split('=')[0].split(' ')[0]
    node_state = line.split('=')[1].split(' ')
    while '' in node_state: node_state.remove('')
    if node_state[0] == '0' or node_state[0] == '1' or node_state[0] == '2':
        node_state = node_state[0]
        key = node + '=' + node_state
        temp = line.split('=')
        temp.pop(0)
        temp2 = ''.join(temp[i] + '=' for i in range(len(temp)))[:-1].replace('if', '').replace('AND', '&&').replace('and', '&&').replace('XOR', '##').replace('OR', '||').replace('or', '||').replace('unless', '!!').replace(' ', '')[1:]
        if 'regardlessof' in temp2:
            temp2 = temp2.split('regardlessof')[0]
        if not temp2 == 'otherwise':
            if '||' in temp2: # OR
                or_node = temp2.split('||')
                other_state = or_node[1][0]
                or_node = or_node[0][len(or_node[0]) - 7:]
                temp2 = temp2.replace(or_node + '||' + other_state, or_node + ' OR ' + or_node[:-1] + other_state)
            if '##' in temp2: # XOR
                xor_node_a = temp2.split('##')
                xor_node_b = xor_node_a[1].split(')')[0]
                xor_node_a = xor_node_a[0][1:]
                temp2 = temp2.replace(xor_node_a + '##' + xor_node_b, '(' + xor_node_a + ' AND NOT ' + xor_node_b + ') OR (NOT ' + xor_node_a + ' AND ' + xor_node_b + ')')
            temp2 = temp2.replace('&&', ' AND ').replace('!!', ' AND NOT ')
            if rdict[key] == '':
                rdict[key] = temp2
            else:
                rdict[key] = '(' + rdict[key] + ') OR (' + temp2 + ')'
    else:
        regulator_node = node_state[0]
        key_prefix = node + '='
        i = 0
        while key_prefix + str(i) in rdict.keys():
            key = key_prefix + str(i)
            if rdict[key] == '':
                rdict[key] = regulator_node + '=' + str(i)
            i += 1

bad_keys = []
for key in rdict.keys():
    if key.split('=')[1] == '0':
        if rdict[key] == '':
            bad_keys.append(key)
        else:
            key_pre = key.split('=')[0] + '='
            i = 1
            while key_pre + str(i) in rdict.keys():
                other_key = key_pre + str(i)
                zero_rule = ' AND NOT (' + rdict[key] + ')'
                rdict[other_key] = '(' + rdict[other_key] + ')' + zero_rule
                i += 1
            bad_keys.append(key)

for key in bad_keys:
    rdict.pop(key)

sort = sorted(rdict)
g = open(pmid + '.txt', 'w')
for key in sort:
    g.write(key + ' :\t' + rdict[key] + '\n')
g.close()