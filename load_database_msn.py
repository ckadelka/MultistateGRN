#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 12:55:34 2025

@author: benco

#Version 3
"""

import math
import numpy as np
import multistate_toolbox as mst

def find_all_indices(arr, el):
    res = []
    for i, a in enumerate(arr):
        if a == el:
            res.append(i)
    if res == []:
        raise ValueError('The element is not in the array')
    return res

def text_to_MSN(folder, file_name, separator = ":", original_equality = "=", original_not = "NOT", original_and = "AND", original_or = "OR", new_equality = " == ", new_not = " not ", new_and = " and ", new_or = " or ", max_degree = 15, ADD_CONSTANTS_TO_NETWORK = True, max_node_count = 10000):
    # Get the text
    f = open(folder + file_name,'r')
    text = f.read()
    text = text.replace('\t',' ').replace('(',' ( ').replace(')',' ) ')
    tvec = text.splitlines()
    f.close()
    
    # remove all empty lines
    while '' in tvec:
        tvec.remove('')
    assert len(tvec) > 0, 'The file is empty'
    
    # the first line is inherently unique
    node_count = 1
    
    # get the count of unique nodes -> non-boolean nodes exist across multiple adjacent lines
    for i, line in enumerate(tvec):
        if i > 0 and not line.split(original_equality)[0] == tvec[i - 1].split(original_equality)[0]:
            node_count += 1
    assert node_count <= max_node_count, 'n=' + str(node_count) + ' > max_n=' + str(max_node_count)
    
    # Determining Variables
    variables = []
    for line in tvec:
        line_var = line[0:line.find(original_equality)].replace(' ', '')
        if not line_var in variables:
            variables.append(line_var)
    
    # Determining Constants
    constants_and_variables = []
    for line in tvec:
        linesplit = line.split(' ')
        for el in linesplit:
            el = el.split(original_equality)[0]
            if el not in ['(', ')', '+', '*', separator, original_not, original_and, original_or, '', ' ']:
                constants_and_variables.append(el)
    constants = list(set(constants_and_variables)-set(variables))
    
    dict_vars_and_consts = dict({original_equality : new_equality, original_not : new_not, original_and : new_and, original_or : new_or})
    dict_vars_and_consts.update(dict(list(zip(variables, ["x[%i]" % i for i in range(len(variables))]))))
    # append constants to the end
    dict_vars_and_consts.update(list(zip(constants, ["x[%i]" % i for i in range(len(variables), len(set(constants_and_variables)))])))
    
    # switch tvec to use the new representation instead of the original ones
    for i, line in enumerate(tvec):
        linesplit = line.split(' ')
        for j, el in enumerate(linesplit):
            if el not in ['(', ')', '+', '*', separator, new_equality.strip(' '), new_not.strip(' '), new_and.strip(' '), new_or.strip(' '), '', ' ']:
                if original_equality in el:
                    equality_split = el.split(original_equality)
                    linesplit[j] = dict_vars_and_consts[equality_split[0]] + new_equality + equality_split[1]
                else:
                    linesplit[j] = dict_vars_and_consts[el]
        tvec[i] = ' '.join(linesplit)
    
    # get the right side of the state equality, the half the describes the rules
    tvec_rule = []
    for i in range(len(tvec)):
        tvec_rule.append(tvec[i][tvec[i].find(separator)+len(separator):])
    
    B = []
    # get the bases for every variable
    
    # We have to check every instance of the variable in the file, as there is potential
    # (albeit rarely does it occur) that there exists some gene that has more states
    # than regulatory rules (ex. a quaternary node that has rules for 1 and 2, but not for 3)
    # If only check the left-side, it becomes much faster, but it will inaccurately read
    # some nodes that have the aforementioned affliction
    for i, v in enumerate(variables):
        node_rep = dict_vars_and_consts[v]
        highest_state = 0
        for line in tvec:
            if node_rep in line:
                state_counts = line.split(node_rep + new_equality)
                for j in state_counts:
                    if j.split(' ')[0].isdigit():
                        state = int(j.split(' ')[0]) + 1
                        if state > highest_state:
                            highest_state = state;
        B.append(highest_state)
    
    # get the bases for every constant
    
    # constants are sporatically placed throughout the file - as they do not have regulatory
    # rules they have no easy way to determine what their highest state count is aside from
    # checking every entry in the file
    for i, v in enumerate(constants):
        node_rep = dict_vars_and_consts[v]
        highest_state = 0
        for line in tvec_rule:
            if node_rep in line:
                state_counts = line.split(node_rep + new_equality)
                for j in state_counts:
                    if j.split(' ')[0].isdigit():
                        state = int(j.split(' ')[0]) + 1
                        if state > highest_state:
                            highest_state = state
        B.append(highest_state)
    
    # convert B to a numpy array from a list
    B = np.array(B, int)
    
    I = []
    # get the regulators for every variable
    for i in range(len(tvec)):
        regulatee = int(tvec[i][find_all_indices(tvec[i], '[')[0] + 1:find_all_indices(tvec[i], ']')[0]])
        idx_open = find_all_indices(tvec_rule[i], '[')
        idx_close = find_all_indices(tvec_rule[i], ']')
        dummy = np.sort(np.array(list(map(int, list(set([tvec_rule[i][(begin+1):end] for begin, end in zip(idx_open, idx_close)]))))))
        if regulatee < len(I):
            I[regulatee] = np.sort(list(set(np.concatenate((I[regulatee], dummy), dtype = int))))
        else:
            I.append(dummy)
    
    degree = list(map(len,I))
    
    # create the right side of the truth table for every variable
    F = []
    for i in range(node_count):
        state_count = math.prod(B[I[i]])
        f = np.zeros(state_count, int)
        F.append(f)
    for i, line in enumerate(tvec):
        node_idx = int(tvec[i][find_all_indices(tvec[i], '[')[0] + 1:find_all_indices(tvec[i], ']')[0]])
        state = int(line.split(str(node_idx) + ']' + new_equality)[1].split(' ')[0])
        x = np.zeros(I[node_idx][len(I[node_idx]) - 1] + 1, int)
        for j, _ in enumerate(F[node_idx]):
            vecj = mst.dec2multi(j, B[I[node_idx]])
            for k, I_v in enumerate(I[node_idx]):
                x[I_v] = vecj[k]
            if eval(tvec_rule[i]):
                F[node_idx][j] = state
    
    if ADD_CONSTANTS_TO_NETWORK:
        # add self-regulation to F and I for every constant
        for i in range(len(constants)):
            f = []
            for b in range(B[i + len(variables)]):
                f.append(b)
            F.append(np.array(f))
            I.append(np.array([len(variables) + i]))
            degree.append(1)
    
    return F, I, B, degree, variables, constants

def create_text_MSN_from_text_BN(folder, file_name, file_name_appendage = '_as_multi', original_separator = '=', original_not = "NOT", original_and = "AND", original_or = "OR", new_separator = ':', equality = "=", new_not = "NOT", new_and = "AND", new_or = "OR"):
    f = open(folder + file_name,'r')
    text = f.read()
    f.close()
    
    text = text.replace('(', ' ( ').replace(')', ' ) ').split('\n')
    mstext = ''
    for i in range(len(text)):
        line = text[i].split(original_separator)
        prefix = line[0].replace(' ', '') + equality + '1 ' + new_separator + '\t'
        rule = line[1].split(' ')
        while '' in rule: rule.remove('')
        for j in range(len(rule)):
            if not (rule[j] == original_and or rule[j] == original_or or rule[j] == original_not or rule[j] == '(' or rule[j] == ')'):
                rule[j] += equality + '1'
        rule = ''.join(rule[j] + ' ' for j in range(len(rule)))
        mstext += (prefix + rule).replace('( ', '(').replace(' )', ')').replace(original_not, new_not).replace(original_and, new_and).replace(original_or, new_or) + '\n'
    
    fullstop_idx = find_all_indices(file_name, '.')[0]
    g = open(folder + file_name[:fullstop_idx] + file_name_appendage + file_name[fullstop_idx:], 'w')
    g.write(mstext)
    g.close()