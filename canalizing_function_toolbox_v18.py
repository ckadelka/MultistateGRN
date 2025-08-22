#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 5 09:53:40 2025

@author: ckadelka

#Version 18

This toolbox allows to 

1) determine if a Boolean function is 
a) constant
b) degenerated
c) canalizing
d) k-canalizing

2) determine the 
a) canalizing depth of a Boolean function 
b) the layer structure of a canaliizing Boolean function

3) generate uniformly at random
a) non-degenerated Boolean functions
a) non-canalizing Boolean functions
c) non-canalizing non-degenerated Boolean functionsr
d) k-canalizing Boolean functions 
e) k-canalizing Boolean functions with a given layer structure
f) Boolean functions with exact canalizing depth k
g) Boolean functions with exact canalizing depth k with a given layer structure

4) generate uniformly at random Boolean networks with specific characterists (in-degree, canalization, strong connectedness)

5) obtain some basic estimates of the magnitude of various subclasses of Boolean functions

6) determine the 
a) absolute bias of a Boolean function
b) average sensitivity of a Boolean function
"""

#13:  proper documentation added, deleted functions that became obsolete
#1.9: new functionality added: calculate feed forward loops and feedback loops
#1.5: added is_kset_canalizing
#1.4: fixed issue with k==0 and EXACT_DEPTH==True in random_BN
#1.3: Python3.7 compatible, kis passed to random_BN can also be [0] for random networks
#1.2: added functionality to randomly create and analyze Boolean networks based on random or canalizing functions
#1.1: fixed a couple issues in is_k_canalizing, is_k_canalizing_return_inputs_outputs_corefunction and get_layerstructure_given_canalizing_outputs_and_corefunction

##Imports

import numpy as np
import matplotlib.pyplot as plt
import itertools
import networkx as nx
import random

# import sympy
# import pandas as pd

from collections import defaultdict
from matplotlib import colors as mcolors
from scipy.special import binom

try:
    import cana.boolean_node
    LOADED_CANA=True
except:
    LOADED_CANA=False

## 0) Basics, helper functions   
import numpy as np
import matplotlib.pyplot as plt
import itertools
import networkx as nx
import random

# import sympy
# import pandas as pd

from collections import defaultdict
from matplotlib import colors as mcolors
from scipy.special import binom

try:
    import cana.boolean_node
    LOADED_CANA = True
except:
    LOADED_CANA = False

## 0) Basics, helper functions   

def bin2dec(binary_vector):
    """
    Convert a binary vector to an integer.

    Parameters:
        binary_vector (list): List containing binary digits (0 or 1).

    Returns:
        int: Integer value converted from the binary vector.
    """
    binary_string = ''.join(str(bit) for bit in binary_vector)
    return int(binary_string, 2)


def dec2bin(integer_value, num_bits):
    """
    Convert an integer to a binary vector.

    Parameters:
        integer_value (int): Integer value to be converted.
        num_bits (int): Number of bits in the binary representation.

    Returns:
        list: List containing binary digits (0 or 1).
    """
    binary_string = bin(integer_value)[2:].zfill(num_bits)
    return [int(bit) for bit in binary_string]


def find_all_indices(array, value):
    """
    Return a list of all indices in the array that equal a given value.

    Parameters:
        array (list): The list in which to search.
        value (any): The value to search for.

    Returns:
        list: List of indices where the value is found.

    Raises:
        ValueError: If the value is not found in the array.
    """
    res = []
    for i, a in enumerate(array):
        if a == value:
            res.append(i)
    if res == []:
        raise ValueError('The value is not in the array at all.')
    return res


def edgelist_to_I(edgelist):
    """
    Convert an edge list describing regulations to an incidence-like list.

    Parameters:
        edgelist (list of lists/tuples): An m x 2 array where each element is a pair 
            [regulator, regulated_node].

    Returns:
        tuple: 
            - I (list of lists): Each element is a list of indices corresponding to regulators 
              for the given node.
            - var (list): List of all unique variables that appear as regulators and/or regulated nodes.
    """
    regulators = np.array(edgelist)[:, 0]
    targets = np.array(edgelist)[:, 1]
    var = list(set(regulators) | set(targets))
    n_var = len(var)
    dict_var = dict(zip(var, range(n_var)))
    I = [[] for _ in range(n_var)]
    for i in range(len(regulators)):
        I[dict_var[targets[i]]].append(dict_var[regulators[i]])
    return I, var


def bool_to_poly(f, left_side_of_truth_table=[], prefix='x', indices=[]):
    """
    Transform a Boolean function from truth table format to polynomial format in non-reduced DNF.

    Parameters:
        f (list): Boolean function as a vector (list of length 2^n, where n is the number of inputs).
        left_side_of_truth_table (list, optional): The left-hand side of the Boolean truth table 
            (a list of tuples of size 2^n x n). If provided, it speeds up computation.
        prefix (str, optional): Prefix for variable names in the polynomial.
        indices (list, optional): List of indices to use for variable naming. If empty or not matching 
            the required number, defaults to list(range(1, n+1)).

    Returns:
        str: A string representing the Boolean function in disjunctive normal form (DNF).
    """
    len_f = len(f)
    n = int(np.log2(len_f))
    if len(indices) == 0 or len(indices) != n:
        indices = list(range(1, n + 1))
    if left_side_of_truth_table == []:  # to reduce run time, this should be calculated once and then passed as argument
        left_side_of_truth_table = list(itertools.product([0, 1], repeat=n))
    num_values = 2 ** n
    text = []
    for i in range(num_values):
        if f[i] == True:
            monomial = '*'.join([(prefix + '%i' % (j)) if entry == 1 else ('(1-' + prefix + '%i)' % (j)) 
                                  for j, entry in zip(indices, left_side_of_truth_table[i])])
            text.append(monomial)
    if text != []:
        return ' + '.join(text)
    else:
        return '0'


def check_if_empty(my_list):
    """
    Check if the provided list or NumPy array is empty.

    Parameters:
        my_list (list or np.ndarray): The list or array to check.

    Returns:
        bool: True if my_list is empty (or has size 0 for a NumPy array), False otherwise.
    """
    if isinstance(my_list, np.ndarray):
        if my_list.size == 0:
            return True
    elif my_list == []:
        return True
    else:
        return False


def eval_expr(expr, x):
    """
    Evaluate a Boolean expression using the shunting-yard algorithm.

    This function parses and evaluates a Boolean expression provided as a string.
    It supports the operators "not", "and", "or" as well as parentheses for grouping.
    The expression must have tokens separated by spaces.

    Parameters:
        expr (str): The Boolean expression in infix notation.
        x (list): A list of Boolean values (0 or 1) corresponding to variables referenced in the expression.
                  Variables in the expression are expected to be in the format "x[i]" where i is an index.

    Returns:
        int: The Boolean result (0 or 1) obtained from evaluating the expression.
    """
    # Helper function to apply the first operator from the operator stack
    def apply_first_op(op_stack, val_stack):
        assert len(op_stack) > 0
        operator = op_stack.pop()
        if operator == "not":
            val_stack.append(int(not val_stack.pop()))
            return

        val1 = val_stack.pop()
        val2 = val_stack.pop()
        outval = apply_operator(operator, val1, val2)
        val_stack.append(outval)

    # Helper function to apply a specific operator to two values
    def apply_operator(operator, val1, val2):
        if operator == "not":
            return int(not val1)
        elif operator == "and":
            return int(val1 and val2)
        elif operator == "or":
            return int(val1 or val2)
        else:
            print("Error: unrecognized operator:", operator)

    operators = {
        "or": 1,
        "and": 2,
        "not": 3,
        "(": 18,
        ")": 18
    }
    
    op_stack = []
    val_stack = []
    for token in expr.split(' '):
        if token == '':
            continue

        if token.isdigit():
            val_stack.append(int(token))
        elif token not in operators:
            # Expected variable format is "x[i]"
            val = x[int(token[2:-1])]
            val_stack.append(val)
        elif token == '(':
            op_stack.append(token)
        elif token == ')':
            while op_stack[-1] != '(':
                apply_first_op(op_stack, val_stack)
            op_stack.pop()  # Remove the '('
        else:
            while (len(op_stack) > 0 and op_stack[-1] != "(" and 
                   operators[op_stack[-1]] >= operators[token]):
                apply_first_op(op_stack, val_stack)
            op_stack.append(token)
        
    while len(op_stack) > 0:
        apply_first_op(op_stack, val_stack)

    return val_stack[0]


def f_from_expression(expr):
    """
    Extract a Boolean function from a string expression.

    The function converts an input expression into its truth table representation.
    The expression can include Boolean operators and comparisons, and the order of variables
    is determined by their first occurrence in the expression.

    Parameters:
        expr (str): A text string containing an evaluable Boolean expression.
            Examples:
                'A AND NOT B'
                'x1 + x2 + x3 > 1'
                '(x1 + x2 + x3) % 2 == 0'

    Returns:
        tuple:
            - f (list): The right-hand side of the Boolean function (truth table) as a list of length 2**n,
              where n is the number of inputs.
            - var (list): A list of variable names (of length n) in the order they were encountered.
    
    Examples:
        >>> f_from_expression('A AND NOT B')
        ([0, 0, 1, 0], ['A', 'B'])
        
        >>> f_from_expression('x1 + x2 + x3 > 1')
        ([0, 0, 0, 1, 0, 1, 1, 1], ['x1', 'x2', 'x3'])
        
        >>> f_from_expression('(x1 + x2 + x3) % 2 == 0')
        ([1, 0, 0, 1, 0, 1, 1, 0], ['x1', 'x2', 'x3'])
    """
    expr = expr.replace('(', ' ( ').replace(')', ' ) ')
    expr_split = expr.split(' ')
    var = []
    dict_var = dict()
    n_var = 0
    for i, el in enumerate(expr_split):
        if el not in ['',' ','(',')','and','or','not','AND','OR','NOT','&','|','~','+','-','*','%','>','>=','==','<=','<'] and not el.isdigit():
            try:
                new_var = dict_var[el]
            except KeyError:
                new_var = 'x[%i]' % n_var
                dict_var.update({el: new_var})
                var.append(el)
                n_var += 1
            expr_split[i] = new_var
        elif el in ['AND','OR','NOT']:
            expr_split[i] = el.lower()
    expr = ' '.join(expr_split)
    f = []
    for x_val in itertools.product([0, 1], repeat=n_var):
        x_val = list(map(bool, x_val))
        f.append(int(eval(expr)))  # x_val is used implicitly in the eval context
    return f, var
            
def enumerate_hypercube_edges(n):
    """
    Enumerate the edges of an n-dimensional hypercube.

    Each vertex of the n-dimensional hypercube is represented as an n-tuple of binary values.
    An edge connects two vertices that differ in exactly one coordinate.
    Duplicate edges (i.e., (v, neighbor) and (neighbor, v)) are avoided.

    Parameters:
        n (int): Dimension of the hypercube.

    Returns:
        list: A list of tuples, where each tuple (v, neighbor) represents an edge between two vertices.
              Each vertex is represented as an n-tuple of 0s and 1s.
    """
    vertices = list(itertools.product([0, 1], repeat=n))
    edges = []
    
    for v in vertices:
        for i in range(n):
            neighbor = list(v)
            neighbor[i] ^= 1  # Flip the i-th bit.
            neighbor = tuple(neighbor)
            if v < neighbor:  # Avoid duplicate edges.
                edges.append((v, neighbor))
    
    return edges

## 1) Methods to analyze Boolean functions
def is_monotonic(F,GET_DETAILS=False):
    n=int(np.log2(len(F)))
    F = np.array(F)
    monotonic = []
    for i in range(n):
        dummy_add=(2**(n-1-i))
        dummy=np.arange(2**n)%(2**(n-i))//dummy_add
        diff = F[dummy==1]-F[dummy==0]
        min_diff = min(diff)
        max_diff = max(diff)
        if min_diff==0 and max_diff==0:
            monotonic.append('not essential')
        elif min_diff==-1 and max_diff==1:
            monotonic.append('not monotonic')
        elif min_diff>=0 and max_diff==1:
            monotonic.append('increasing')            
        elif min_diff==-1 and max_diff<=0:
            monotonic.append('decreasing')   
    if GET_DETAILS:
        return ('not essential' not in monotonic,monotonic)
    else:
        return 'not essential' not in monotonic


def is_degenerated(f):
    """
    Determine if a Boolean function contains non-essential variables.

    A variable is non-essential if the function's output does not depend on it.

    Parameters:
        f (list): Boolean function represented as a list of length 2^n (truth table), where n is the number of inputs.

    Returns:
        bool: True if f contains at least one non-essential variable, False if all variables are essential.
    """
    len_f = len(f)
    n = int(np.log2(len_f))
    for i in range(n):
        dummy_add = (2**(n-1-i))
        dummy = np.arange(2**n) % (2**(n-i)) // dummy_add
        depends_on_i = False
        for j in range(2**n):
            if dummy[j] == 1:
                continue
            else:
                if f[j] != f[j + dummy_add]:
                    depends_on_i = True
                    break
        if depends_on_i == False:
            return True
    return 
def get_essential_variables(f):
    """
    Determine the indices of essential variables in a Boolean function.

    A variable is essential if changing its value (while holding the others constant) can change the output of f.

    Parameters:
        f (list): Boolean function as a list of length 2^n (truth table), where n is the number of inputs.

    Returns:
        list: List of indices corresponding to the essential variables.
    """
    if len(f) == 0:
        return []
    len_f = len(f)
    n = int(np.log2(len_f))
    essential_variables = list(range(n))
    for i in range(n):
        dummy_add = (2**(n-1-i))
        dummy = np.arange(2**n) % (2**(n-i)) // dummy_add
        depends_on_i = False
        for j in range(2**n):
            if dummy[j] == 1:
                continue
            else:
                if f[j] != f[j + dummy_add]:
                    depends_on_i = True
                    break
        if depends_on_i == False:
            essential_variables.remove(i)
    return essential_variables 


def get_number_essential_variables(f):
    """
    Count the number of essential variables in a Boolean function.

    Parameters:
        f (list): Boolean function as a list of length 2^n (truth table), where n is the number of inputs.

    Returns:
        int: The number of essential variables.
    """
    return len(get_essential_variables(f))


def is_constant(f):
    """
    Check whether a Boolean function is constant.

    Parameters:
        f (list): Boolean function as a list of length 2^n (truth table), where n is the number of inputs.

    Returns:
        bool: True if f is constant (all outputs are 0 or all are 1), False otherwise.
    """
    return sum(f) in [0, len(f)]


def get_symmetry_groups(f, left_side_of_truth_table=[]):
    """
    Determine all symmetry groups of input variables for a Boolean function.

    Two variables are in the same symmetry group if swapping their values does not change the output
    of the function for any input of the other variables.

    Parameters:
        f (list or np.array): Boolean function of length 2^n (truth table), where n is the number of inputs.
        left_side_of_truth_table (optional, array-like): Precomputed left-hand side of the truth table (2^n x n).
            If not provided or if its shape does not match, it will be computed.

    Returns:
        list: A list of lists where each inner list contains indices of variables that form a symmetry group.
    """
    len_f = len(f)
    n = int(np.log2(len_f))
    if left_side_of_truth_table == [] or left_side_of_truth_table.shape[0] != len_f:
        left_side_of_truth_table = np.array(list(itertools.product([0, 1], repeat=n)))
    symmetry_groups = []
    left_to_check = np.ones(n)
    for i in range(n):
        if left_to_check[i] == 0:
            continue
        else:
            symmetry_groups.append([i])
            left_to_check[i] = 0
        for j in range(i + 1, n):
            diff = sum(2**np.arange(n - i - 2, n - j - 2, -1))
            for ii, x in enumerate(left_side_of_truth_table):
                if x[i] != x[j] and x[i] == 0 and f[ii] != f[ii + diff]:
                    break
            else:
                left_to_check[j] = 0
                symmetry_groups[-1].append(j)
    return symmetry_groups

def get_average_sensitivity(f, nsim=10000, EXACT=False, NORMALIZED=True):
    """
    Compute the average sensitivity of a Boolean function.

    The average sensitivity is equivalent to the Derrida value D(F,1) when the update rule is sampled
    from the same space. This function can compute the exact sensitivity by exhaustively iterating over all inputs (if EXACT is True)
    or estimate it via Monte Carlo sampling (if EXACT is False). The result can be normalized by the number of inputs.

    Parameters:
        f (list or np.array): Boolean function (truth table) of length 2^n, where n is the number of inputs.
        nsim (int, optional): Number of random samples (default is 10000, used when EXACT is False).
        EXACT (bool, optional): If True, compute the exact sensitivity by iterating over all inputs; otherwise, use sampling (default).
        NORMALIZED (bool, optional): If True, return the normalized sensitivity (divided by the number of function inputs); otherwise, return the total count.

    Returns:
        float: The (normalized) average sensitivity of the Boolean function.
    """
    if type(f) == list:
        f = np.array(f)
    n = int(np.log2(len(f)))
    num_values = 2**n
    s = 0
    if EXACT:
        left_side_of_truth_table = list(map(np.array, list(itertools.product([0, 1], repeat=n))))
        for ii, X in enumerate(left_side_of_truth_table):
            for i in range(n):
                Y = X.copy()
                Y[i] = 1 - X[i]
                Ydec = bin2dec(Y)
                s += int(f[ii] != f[Ydec])
        if NORMALIZED:
            return s / (num_values * n)
        else:
            return s / num_values
    else:
        for i in range(nsim):
            xdec = np.random.randint(num_values)
            Y = dec2bin(xdec, n)
            index = np.random.randint(n)
            Y[index] = 1 - Y[index]
            Ybin = bin2dec(Y)
            s += int(f[xdec] != f[Ybin])
        if NORMALIZED:
            return s / nsim
        else:
            return n * s / nsim


def get_absolute_bias(f, n=None):
    """
    Compute the absolute bias of a Boolean function.

    The absolute bias is defined as |(sum(f) / 2^(n-1)) - 1|, which quantifies how far the function's output distribution
    deviates from being balanced.

    Parameters:
        f (list or np.array): Boolean function (truth table) of length 2^n.
        n (int, optional): Number of inputs. If not provided, n is inferred from the length of f.

    Returns:
        float: The absolute bias of the Boolean function.
    """
    if n is None:
        n = int(np.log2(len(f)))
    return abs(sum(f) * 1.0 / 2**(n - 1) - 1)



def is_canalizing(f, n=-1):
    """
    Determine if a Boolean function is canalizing.

    A Boolean function f(x_1, ..., x_n) is canalizing if there exists at least one variable x_i and a value a ∈ {0, 1} 
    such that f(x_1, ..., x_i = a, ..., x_n) is constant.

    Parameters:
        f (list or np.array): Boolean function of length 2^n (truth table), where n is the number of inputs.
        n (int, optional): Number of variables. If n == -1, it is inferred from the length of f.

    Returns:
        bool: True if f is canalizing, False otherwise.
    """
    if type(f) == list:
        f = np.array(f)
    if n == -1:
        n = int(np.log2(len(f)))
    desired_value = 2**(n - 1)
    T = np.array(list(itertools.product([0, 1], repeat=n))).T
    A = np.r_[T, 1 - T]
    Atimesf = np.dot(A, f)
    if np.any(Atimesf == desired_value):
        return True
    elif np.any(Atimesf == 0):
        return True
    else:
        return False


def is_kset_canalizing(f, k, n=-1):
    """
    Determine if a Boolean function is k-set canalizing.

    A Boolean function is k-set canalizing if there exists a set of k variables such that setting these variables to specific values
    forces the output of the function, irrespective of the other n - k inputs.

    Parameters:
        f (list or np.array): Boolean function of length 2^n (truth table), where n is the number of inputs.
        k (int): The size of the variable set (with 0 ≤ k ≤ n).
        n (int, optional): Number of variables. If n == -1, it is inferred from the length of f.

    Returns:
        bool: True if f is k-set canalizing, False otherwise.

    References:
        Kadelka, C., Keilty, B., & Laubenbacher, R. (2023). Collectively canalizing Boolean functions.
        Advances in Applied Mathematics, 145, 102475.
    """
    if type(f) == list:
        f = np.array(f)
    if k == 0:
        return is_constant(f)
    if n == -1:
        n = int(np.log2(len(f)))
    desired_value = 2**(n - k)
    T = np.array(list(itertools.product([0, 1], repeat=n))).T
    A = np.r_[T, 1 - T]
    Ak = []
    for i in range(2 * n):
        for j in range(i + 1, 2 * n):
            if j - i == n:
                continue
            else:
                Ak.append(np.bitwise_and(A[i, :], A[j, :]))
    Ak = []
    for indices in itertools.combinations(range(2 * n), k):
        dummy = np.sum(A[np.array(indices), :], 0) == k
        if sum(dummy) == desired_value:
            Ak.append(dummy)
    Ak = np.array(Ak)
    AktimesF = np.dot(Ak, f)
    is_there_canalization = 0 in AktimesF or desired_value in AktimesF
    return is_there_canalization


def get_proportion_of_collectively_canalizing_input_sets(f, k, n=-1, left_side_of_truth_table=[], verbose=False):
    """
    Compute the proportion of k-set canalizing input sets for a Boolean function.

    For a given k, this function calculates the probability that a randomly chosen set of k inputs canalizes the function,
    i.e., forces the output regardless of the remaining variables.

    Parameters:
        f (list or np.array): Boolean function of length 2^n (truth table), where n is the number of inputs.
        k (int): The size of the variable set (0 ≤ k ≤ n).
        n (int, optional): Number of variables. If n == -1, it is inferred from the length of f.
        left_side_of_truth_table (optional, array-like): Precomputed left-hand side of the truth table (2^n x n).
            If not provided, it is computed.
        verbose (bool, optional): If True, prints detailed information about canalizing k-sets.

    Returns:
        float: The proportion of k-set canalizing input sets.
    
    References:
        Kadelka, C., Keilty, B., & Laubenbacher, R. (2023). Collectively canalizing Boolean functions.
        Advances in Applied Mathematics, 145, 102475.
    """
    if type(f) == list:
        f = np.array(f)
    if k == 0:
        return float(is_constant(f))
    if n == -1:
        n = int(np.log2(len(f)))
    desired_value = 2**(n - k)
    if left_side_of_truth_table == []:
        T = np.array(list(itertools.product([0, 1], repeat=n))).T
    else:
        T = np.array(left_side_of_truth_table).T
    Tk = list(itertools.product([0, 1], repeat=k))
    A = np.r_[T, 1 - T]
    Ak = []
    for indices in itertools.combinations(range(n), k):
        for canalizing_inputs in Tk:
            indices_values = np.array(indices) + n * np.array(canalizing_inputs)
            dummy = np.sum(A[indices_values, :], 0) == k
            if sum(dummy) == desired_value:
                Ak.append(dummy)
                if verbose and np.dot(dummy, f) in [0, desired_value]:
                    print(indices, canalizing_inputs, indices_values, np.dot(dummy, f))
            elif verbose:
                print(indices, canalizing_inputs, sum(dummy), 'a')
    Ak = np.array(Ak)
    is_there_canalization = np.in1d(np.dot(Ak, f), [0, desired_value])
    return sum(is_there_canalization) / len(is_there_canalization)


def get_canalizing_strength(f, left_side_of_truth_table=[]):
    """
    Compute the canalizing strength of a Boolean function via exhaustive enumeration.

    The canalizing strength is defined as a weighted average of the proportions of k-set canalizing inputs for k = 1 to n-1.
    It is 0 for minimally canalizing functions (e.g., Boolean parity functions) and 1 for maximally canalizing functions
    (e.g., nested canalizing functions with one layer).

    Parameters:
        f (list or np.array): Boolean function of length 2^n (truth table), where n is the number of inputs.
        left_side_of_truth_table (optional, array-like): Precomputed left-hand side of the truth table (2^n x n).

    Returns:
        tuple:
            - float: The canalizing strength of f.
            - list: A list of the k-set canalizing proportions for k = 1, 2, ..., n-1.
    
    References:
        Kadelka, C., Keilty, B., & Laubenbacher, R. (2023). Collectively canalizing Boolean functions.
        Advances in Applied Mathematics, 145, 102475.
    """
    nfloat = np.log2(len(f))
    n = int(nfloat)
    assert abs(n - nfloat) < 1e-10, "f needs to be of length 2^n for some n > 1"
    assert n > 1, "Canalizing strength is only defined for Boolean functions with n > 1 inputs"
    res = []
    for k in range(1, n):
        res.append(get_proportion_of_collectively_canalizing_input_sets(f, k, n, left_side_of_truth_table=left_side_of_truth_table))
    return np.mean(np.multiply(res, 2**np.arange(1, n) / (2**np.arange(1, n) - 1))), res


def compute_exact_kset_canalizing_proportion_for_ncf_with_specific_layerstructure(k, layerstructure_NCF):
    """
    Compute the exact k-set canalizing proportion for a nested canalizing function (NCF) with a specific layer structure.

    This function implements Theorem 3.3 from [1] and computes the exact proportion of k-set canalizing inputs for an NCF
    characterized by its layer structure.

    Parameters:
        k (int): The size of the variable set (0 ≤ k ≤ n) for which the canalizing proportion is computed.
        layerstructure_NCF (list): List of integers [k_1, ..., k_r] describing the number of variables in each layer of an NCF.
            Each k_i must be at least 1, and the last layer must have at least 2 variables unless n == 1.

    Returns:
        float: The exact k-set canalizing proportion for the NCF with the provided layer structure.
    
    References:
        [1] Kadelka, C., Keilty, B., & Laubenbacher, R. (2023). Collectively canalizing Boolean functions.
            Advances in Applied Mathematics, 145, 102475.
    """
    r = len(layerstructure_NCF)
    n = sum(layerstructure_NCF)
    assert min(layerstructure_NCF) >= 1 and (layerstructure_NCF[-1] >= 2 or n == 1), \
           "Each layer must contain at least one variable (the last layer at least two unless n == 1)"
    magnitudes = []
    for t in range(r):
        number_of_input_sets = 0
        for c in range(1, min(k - sum(layerstructure_NCF[:t][::-2]), layerstructure_NCF[t]) + 1):
            for d in range(0, min(k - sum(layerstructure_NCF[:t][::-2]) - c, sum(layerstructure_NCF[:max(0, t - 1)][::-2])) + 1):
                binom1 = binom(layerstructure_NCF[t], c)
                binom2 = binom(sum(layerstructure_NCF[:max(0, t - 1)][::-2]), d)
                binom3 = binom(n - sum(layerstructure_NCF[:t + 1]), k - sum(layerstructure_NCF[:t][::-2]) - c - d)
                number_of_inputs_that_canalize_for_selected_variable_set = sum([2**(k - sum(layerstructure_NCF[:t][::-2]) - j - d) for j in range(1, c + 1)])
                number_of_input_sets += binom1 * binom2 * binom3 * number_of_inputs_that_canalize_for_selected_variable_set
        magnitudes.append(number_of_input_sets)
    # For the case where the non-canalizing output value can be reached in the evaluation process, add:
    if k >= sum(layerstructure_NCF[-1::-2]):
        magnitudes.append(binom(n - sum(layerstructure_NCF[-1::-2]), k - sum(layerstructure_NCF[-1::-2])))
    else:
        magnitudes.append(0)
    return sum(magnitudes) / (2**k * binom(n, k))


def is_k_canalizing(f, k, n=-1):
    """
    Determine if a Boolean function is k-canalizing.

    A Boolean function is k-canalizing if it has at least k conditionally canalizing variables.
    This is checked recursively: after fixing a canalizing variable (with a fixed canalizing input that forces the output),
    the subfunction (core function) must itself be canalizing for the next variable, and so on.

    Parameters:
        f (list or np.array): Boolean function of length 2^n (truth table), where n is the number of inputs.
        k (int): The desired canalizing depth (0 ≤ k ≤ n). Note: every function is 0-canalizing.
        n (int, optional): Number of variables. If n == -1, it is inferred from the length of f.

    Returns:
        bool: True if f is k-canalizing, False otherwise.
    
    References:
        He, Q., & Macauley, M. (2016). Stratification and enumeration of Boolean functions by canalizing depth.
            Physica D: Nonlinear Phenomena, 314, 1-8.
        Dimitrova, E., Stigler, B., Kadelka, C., & Murrugarra, D. (2022). Revealing the canalizing structure of Boolean functions:
            Algorithms and applications. Automatica, 146, 110630.
    """
    if k > n:
        return False
    if k == 0:
        return True
    if n == -1:
        n = int(np.log2(len(f)))
    w = sum(f)  # Hamming weight of f
    if w == 0 or w == 2**n:  # constant function
        return False
    if type(f) == list:
        f = np.array(f)
    desired_value = 2**(n - 1)
    T = np.array(list(itertools.product([0, 1], repeat=n))).T
    A = np.r_[T, 1 - T]
    try:  # check for canalizing output 1
        index = list(np.dot(A, f)).index(desired_value)
        new_f = f[np.where(A[index] == 0)[0]]
        return is_k_canalizing(new_f, k - 1, n - 1)
    except ValueError:
        try:  # check for canalizing output 0
            index = list(np.dot(A, 1 - f)).index(desired_value)
            new_f = f[np.where(A[index] == 0)[0]]
            return is_k_canalizing(new_f, k - 1, n - 1)
        except ValueError:
            return False


def is_k_canalizing_return_inputs_outputs_corefunction(f, k, n, can_inputs=np.array([], dtype=int), can_outputs=np.array([], dtype=int)):
    """
    Determine if a Boolean function is k-canalizing and return associated canalizing data.

    This function recursively checks whether f is k-canalizing and returns:
      - A boolean indicating success.
      - The canalizing input values.
      - The canalized output values.
      - The core function that remains after removing the canalizing variables.

    Parameters:
        f (list or np.array): Boolean function of length 2^n (truth table), where n is the number of inputs.
        k (int): The canalizing depth to check.
        n (int): Number of variables.
        can_inputs (np.array, optional): Accumulated canalizing input values (default is an empty array).
        can_outputs (np.array, optional): Accumulated canalized output values (default is an empty array).

    Returns:
        tuple: A tuple containing:
            - bool: True if f is k-canalizing, False otherwise.
            - np.array: Array of canalizing input values.
            - np.array: Array of canalized output values.
            - np.array: The core function (remaining truth table) after canalizing variables are removed.
    
    References:
        He, Q., & Macauley, M. (2016). Stratification and enumeration of Boolean functions by canalizing depth.
            Physica D: Nonlinear Phenomena, 314, 1-8.
        Dimitrova, E., Stigler, B., Kadelka, C., & Murrugarra, D. (2022). Revealing the canalizing structure of Boolean functions:
            Algorithms and applications. Automatica, 146, 110630.
    """
    if k == 0:
        return (True, can_inputs, can_outputs, f)
    w = sum(f)
    if w == 0 or w == 2**n:  # constant function
        return (False, can_inputs, can_outputs, f)
    if type(f) == list:
        f = np.array(f)
    desired_value = 2**(n - 1)
    T = np.array(list(itertools.product([0, 1], repeat=n))).T
    A = np.r_[T, 1 - T]
    try:  # check for canalizing output 1
        index = list(np.dot(A, f)).index(desired_value)
        new_f = f[np.where(A[index] == 0)[0]]
        return is_k_canalizing_return_inputs_outputs_corefunction(new_f, k - 1, n - 1,
                                                                  np.append(can_inputs, int(index < n)),
                                                                  np.append(can_outputs, 1))
    except ValueError:
        try:  # check for canalizing output 0
            index = list(np.dot(A, 1 - f)).index(desired_value)
            new_f = f[np.where(A[index] == 0)[0]]
            return is_k_canalizing_return_inputs_outputs_corefunction(new_f, k - 1, n - 1,
                                                                      np.append(can_inputs, int(index < n)),
                                                                      np.append(can_outputs, 0))
        except ValueError:
            return (False, can_inputs, can_outputs, f)


def is_k_canalizing_return_inputs_outputs_corefunction_order(f, k, n, can_inputs=np.array([], dtype=int),
                                                            can_outputs=np.array([], dtype=int), can_order=np.array([], dtype=int),
                                                            variables=[]):
    """
    Determine if a Boolean function is k-canalizing and return canalizing data including variable order.

    This function extends the k-canalizing check by additionally returning the order (indices) of the canalizing variables.
    It recursively collects:
      - Canalizing input values.
      - Canalized output values.
      - The core function after removing the canalizing layers.
      - The order of the canalizing variables.

    Parameters:
        f (list or np.array): Boolean function of length 2^n (truth table), where n is the number of inputs.
        k (int): The canalizing depth to check.
        n (int): Number of variables.
        can_inputs (np.array, optional): Accumulated canalizing input values.
        can_outputs (np.array, optional): Accumulated canalized output values.
        can_order (np.array, optional): Accumulated order (indices) of canalizing variables.
        variables (list, optional): List of variable indices. If empty, defaults to range(n).

    Returns:
        tuple: A tuple containing:
            - bool: True if f is k-canalizing, False otherwise.
            - np.array: Array of canalizing input values.
            - np.array: Array of canalized output values.
            - np.array: The core function (remaining truth table) after removing canalizing variables.
            - np.array: Array of indices indicating the order of canalizing variables.
    
    References:
        He, Q., & Macauley, M. (2016). Stratification and enumeration of Boolean functions by canalizing depth.
            Physica D: Nonlinear Phenomena, 314, 1-8.
        Dimitrova, E., Stigler, B., Kadelka, C., & Murrugarra, D. (2022). Revealing the canalizing structure of Boolean functions:
            Algorithms and applications. Automatica, 146, 110630.
    """
    if k == 0:
        return (True, can_inputs, can_outputs, f, can_order)
    w = sum(f)
    if w == 0 or w == 2**n:  # constant function
        return (False, can_inputs, can_outputs, f, can_order)
    if type(variables) == np.ndarray:
        variables = list(variables)
    if variables == []:
        variables = list(range(n))
    if type(f) == list:
        f = np.array(f)
    desired_value = 2**(n - 1)
    T = np.array(list(itertools.product([0, 1], repeat=n))).T
    A = np.r_[T, 1 - T]
    try:  # check for canalizing output 0
        index = list(np.dot(A, 1 - f)).index(desired_value)
        newF = f[np.where(A[index] == 0)[0]]
        variable = variables.pop(index % n)
        return is_k_canalizing_return_inputs_outputs_corefunction_order(newF, k - 1, n - 1,
                                                                        np.append(can_inputs, int(index < n)),
                                                                        np.append(can_outputs, 0),
                                                                        np.append(can_order, variable),
                                                                        variables)
    except ValueError:
        try:  # check for canalizing output 1
            index = list(np.dot(A, f)).index(desired_value)
            newF = f[np.where(A[index] == 0)[0]]
            variable = variables.pop(index % n)
            return is_k_canalizing_return_inputs_outputs_corefunction_order(newF, k - 1, n - 1,
                                                                            np.append(can_inputs, int(index < n)),
                                                                            np.append(can_outputs, 1),
                                                                            np.append(can_order, variable),
                                                                            variables)
        except ValueError:
            return (False, can_inputs, can_outputs, f, can_order)


def find_layers(f, can_inputs=np.array([], dtype=int), can_outputs=np.array([], dtype=int),
                can_order=np.array([], dtype=int), variables=[], depth=0, number_layers=0):
    """
    Determine the canalizing layer structure of a Boolean function.

    This function decomposes a Boolean function into its canalizing layers (standard monomial form)
    by recursively identifying and removing conditionally canalizing variables.
    The output includes the canalizing depth, the number of layers, the canalizing inputs and outputs,
    the core polynomial, and the order of the canalizing variables.

    Parameters:
        f (list or np.array): Boolean function of length 2^n (truth table), where n is the number of inputs.
        can_inputs (np.array, optional): Accumulated canalizing input values (for recursion).
        can_outputs (np.array, optional): Accumulated canalized output values (for recursion).
        can_order (np.array, optional): Accumulated indices of canalizing variables (for recursion).
        variables (list, optional): List of variable indices. If empty, defaults to range(n).
        depth (int, optional): Current canalizing depth (for recursion); default is 0.
        number_layers (int, optional): Current number of layers identified (for recursion); default is 0.

    Returns:
        tuple: A tuple containing:
            - int: Canalizing depth (number of conditionally canalizing variables).
            - int: Number of distinct canalizing layers.
            - np.array: Array of canalizing input values.
            - np.array: Array of canalized output values.
            - np.array: The core polynomial (truth table) after removing canalizing variables.
            - np.array: Array of indices representing the order of canalizing variables.
    
    References:
        He, Q., & Macauley, M. (2016). Stratification and enumeration of Boolean functions by canalizing depth.
            Physica D: Nonlinear Phenomena, 314, 1-8.
        Dimitrova, E., Stigler, B., Kadelka, C., & Murrugarra, D. (2022). Revealing the canalizing structure of Boolean functions:
            Algorithms and applications. Automatica, 146, 110630.
    """
    n = int(np.log2(len(f)))
    w = sum(f)
    if w == 0 or w == 2**n:  # constant function
        return (depth, number_layers, can_inputs, can_outputs, f, can_order)
    if type(variables) == np.ndarray:
        variables = list(variables)
    if variables == []:
        variables = list(range(n))
    if type(f) == list:
        f = np.array(f)
    desired_value = 2**(n - 1)
    T = np.array(list(itertools.product([0, 1], repeat=n))).T
    A = np.r_[T, 1 - T]

    indices1 = np.where(np.dot(A, f) == desired_value)[0]
    indices0 = np.where(np.dot(A, 1 - f) == desired_value)[0]
    if len(indices1) > 0:
        sorted_order = sorted(range(len(indices1)), key=lambda x: (indices1 % n)[x])
        inputs = (1 - indices1 // n)[np.array(sorted_order)]
        outputs = np.ones(len(indices1), dtype=int)
        new_canalizing_variables = []
        for index in np.sort(indices1 % n)[::-1]:
            new_canalizing_variables.append(variables.pop(index))
        new_canalizing_variables.reverse()
        new_f = f[np.sort(list(set.intersection(*[] + [set(np.where(A[index] == 0)[0]) for index, INPUT in zip(indices1, inputs)])))]
        return find_layers(new_f, np.append(can_inputs, inputs), np.append(can_outputs, outputs),
                           np.append(can_order, new_canalizing_variables), variables, depth + len(new_canalizing_variables),
                           number_layers + 1)
    elif len(indices0):
        sorted_order = sorted(range(len(indices0)), key=lambda x: (indices0 % n)[x])
        inputs = (1 - indices0 // n)[np.array(sorted_order)]
        outputs = np.zeros(len(indices0), dtype=int)
        new_canalizing_variables = []
        for index in np.sort(indices0 % n)[::-1]:
            new_canalizing_variables.append(variables.pop(index))
        new_canalizing_variables.reverse()
        new_f = f[np.sort(list(set.intersection(*[] + [set(np.where(A[index] == 0)[0]) for index, INPUT in zip(indices0, inputs)])))]
        return find_layers(new_f, np.append(can_inputs, inputs), np.append(can_outputs, outputs),
                           np.append(can_order, new_canalizing_variables), variables, depth + len(new_canalizing_variables),
                           number_layers + 1)
    else:
        return (depth, number_layers, can_inputs, can_outputs, f, can_order)

## 2) Put everything together to obtain canalizing depth, layer structure, canalized outputs, canalizing inputs as well as core function (could also calculate canalizing variables in future versions but I don't see a need)
if LOADED_CANA:
    def get_input_redundancy(f, n=-1):
        """
        Compute the input redundancy of a Boolean function.

        The input redundancy quantifies how many inputs are not required to determine the function’s output.
        Constant functions have an input redundancy of 1 (none of the inputs are needed), whereas parity functions have an input redundancy of 0 (all inputs are necessary).

        Parameters:
            f (list or np.array): Boolean function of length 2^n (truth table), where n is the number of inputs.
            n (int, optional): Number of inputs. If n == -1, it is inferred from len(f).

        Returns:
            float: Normalized input redundancy in the interval [0, 1].

        References:
            [1] Marques-Pita, M., & Rocha, L. M. (2013). Canalization and control in automata networks: body segmentation in Drosophila melanogaster. PloS One, 8(3), e55946.
            [2] Correia, R. B., Gates, A. J., Wang, X., & Rocha, L. M. (2018). CANA: a python package for quantifying control and canalization in Boolean networks. Frontiers in Physiology, 9, 1046.
        """
        if n == -1:
            n = int(np.log2(len(f)))
        return cana.boolean_node.BooleanNode(k=n, outputs=f).input_redundancy()

    def get_edge_effectiveness(f, n=-1):
        """
        Compute the edge effectiveness for each regulator of a Boolean function.

        Edge effectiveness measures how much flipping a given input (regulator) influences the output.
        Non-essential inputs have an effectiveness of 0, whereas inputs that always flip the output when toggled have an effectiveness of 1.

        Parameters:
            f (list or np.array): Boolean function of length 2^n (truth table), where n is the number of inputs.
            n (int, optional): Number of inputs. If n == -1, it is inferred from len(f).

        Returns:
            list: A list of n floats in [0, 1] representing the edge effectiveness for each input.

        References:
            [1] Marques-Pita, M., & Rocha, L. M. (2013). Canalization and control in automata networks: body segmentation in Drosophila melanogaster. PloS One, 8(3), e55946.
            [2] Correia, R. B., Gates, A. J., Wang, X., & Rocha, L. M. (2018). CANA: a python package for quantifying control and canalization in Boolean networks. Frontiers in Physiology, 9, 1046.
        """
        if n == -1:
            n = int(np.log2(len(f)))
        return cana.boolean_node.BooleanNode(k=n, outputs=f).edge_effectiveness()


def get_canalizing_depth_inputs_outputs_corefunction(f):
    """
    (Obsolete) Retrieve the canalizing depth, canalizing inputs, canalized outputs, and core function of a Boolean function.

    This function is maintained for backward compatibility. It is recommended to use find_layers(f) instead.
    The canalizing depth is determined by recursively extracting canalizing variables from the function.

    Parameters:
        f (list or np.array): Boolean function of length 2^n (truth table), where n is the number of inputs.

    Returns:
        tuple: A tuple (n, depth, can_inputs, can_outputs, corefunction), where:
            - n (int): The number of inputs.
            - depth (int): The canalizing depth (number of canalizing inputs).
            - can_inputs (list): List of canalizing input values.
            - can_outputs (list): List of canalized output values.
            - corefunction (list or np.array): The remaining core function after removing canalizing variables.
    """
    n = int(np.log2(len(f)))
    (NESTED_CANALIZING, can_inputs, can_outputs, corefunction) = is_k_canalizing_return_inputs_outputs_corefunction(f, n, n)
    return (n, len(can_inputs), can_inputs, can_outputs, corefunction)


def get_canalizing_depth_inputs_outputs_corefunction_order(f, variables=[]):
    """
    (Obsolete) Retrieve the canalizing depth, canalizing inputs, canalized outputs, core function, and the order of canalizing variables.

    This function is maintained for backward compatibility. It is recommended to use find_layers(f) instead.
    The canalizing depth and related information are computed by recursively extracting canalizing variables.

    Parameters:
        f (list or np.array): Boolean function of length 2^n (truth table), where n is the number of inputs.
        variables (list, optional): List of variable indices. If not provided, it defaults to range(n).

    Returns:
        tuple: A tuple (n, depth, can_inputs, can_outputs, corefunction, can_order), where:
            - n (int): The number of inputs.
            - depth (int): The canalizing depth (number of canalizing inputs).
            - can_inputs (list): List of canalizing input values.
            - can_outputs (list): List of canalized output values.
            - corefunction (list or np.array): The remaining core function after removing canalizing variables.
            - can_order (list): The order of canalizing variable indices.
    """
    n = int(np.log2(len(f)))
    (NESTED_CANALIZING, can_inputs, can_outputs, corefunction, can_order) = \
        is_k_canalizing_return_inputs_outputs_corefunction_order(f, n, n, variables=variables)
    return (n, len(can_inputs), can_inputs, can_outputs, corefunction, can_order)


def get_layerstructure_given_canalizing_outputs_and_corefunction(can_outputs, core_polynomial, n=-1):
    """
    Compute the canalizing layer structure of a Boolean function given its canalized outputs and core polynomial.

    Two consecutive canalizing variables belong to the same layer if they have the same canalized output, and to different layers otherwise.
    The resulting layer structure is a list [k_1, ..., k_r] indicating the number of variables in each canalizing layer.
    For nested canalizing functions (NCFs) with n > 1, the last layer must have at least two variables.

    Parameters:
        can_outputs (list): List of all canalized output values of the function.
        core_polynomial (list or np.array): Core function (or polynomial) of length 2^(n - depth) after removing canalizing variables.
        n (int, optional): Number of inputs of the original function. If n == -1, it is inferred as int(np.log2(len(core_polynomial))) + depth.

    Returns:
        list: A list [k_1, ..., k_r] describing the number of variables in each canalizing layer.
              Each k_i ≥ 1, and if the function is an NCF (i.e., sum(k_i) == n), then the last layer k_r ≥ 2 (unless n == 1).

    References:
        [1] Kadelka, C., Kuipers, J., & Laubenbacher, R. (2017). The influence of canalization on the robustness of Boolean networks.
            Physica D: Nonlinear Phenomena, 353, 39-47.
        [2] Dimitrova, E., Stigler, B., Kadelka, C., & Murrugarra, D. (2022). Revealing the canalizing structure of Boolean functions:
            Algorithms and applications. Automatica, 146, 110630.
    """
    depth = len(can_outputs)
    if depth == 0:
        return []
    if n == -1:
        n = int(np.log2(len(core_polynomial))) + depth
    assert depth != n - 1, ("len(can_outputs) == n-1, which is impossible because the last variable "
                            "must also be canalizing in this case.")
    if depth == n and n > 1:  # For Boolean NCFs, the last layer must have at least two variables.
        can_outputs[-1] = can_outputs[-2]
    elif is_constant(core_polynomial) and depth > 1:  # Exceptional case: last layer needs to be size ≥ 2.
        can_outputs[-1] = can_outputs[-2]
    layerstructure = []
    size_of_layer = 1
    for i in range(1, depth):
        if can_outputs[i] == can_outputs[i - 1]:
            size_of_layer += 1
        else:
            layerstructure.append(size_of_layer)
            size_of_layer = 1
    layerstructure.append(size_of_layer)
    return layerstructure


def get_layerstructure_of_an_NCF_given_its_Hamming_weight(n, w):
    """
    Compute the canalizing layer structure of a nested canalizing function (NCF) given its Hamming weight.

    There exists a bijection between the Hamming weight (with w equivalent to 2^n - w) and the canalizing layer structure of an NCF.
    The layer structure is represented as [k_1, ..., k_r], where each k_i ≥ 1 and, for n > 1, the last layer k_r ≥ 2.

    Parameters:
        n (int): Number of inputs (variables) of the NCF.
        w (int): Odd Hamming weight of the NCF, i.e., the number of 1s in the 2^n-vector representation of the function.

    Returns:
        tuple: A tuple (r, layerstructure_NCF), where:
            - r (int): The number of canalizing layers.
            - layerstructure_NCF (list): A list [k_1, ..., k_r] describing the number of variables in each layer.

    References:
        Kadelka, C., Kuipers, J., & Laubenbacher, R. (2017). The influence of canalization on the robustness of Boolean networks.
        Physica D: Nonlinear Phenomena, 353, 39-47.
    """
    if w == 1:
        r = 1
        layerstructure_NCF = [n]
    else:
        assert type(w) == int or type(w) == np.int64, 'Hamming weight must be an integer'
        assert 1 <= w <= 2**n - 1, 'Hamming weight w must satisfy 1 <= w <= 2^n - 1'
        assert w % 2 == 1, 'Hamming weight must be an odd integer since all NCFs have an odd Hamming weight.'
        w_bin = dec2bin(w, n)
        current_el = w_bin[0]
        layerstructure_NCF = [1]
        for el in w_bin[1:-1]:
            if el == current_el:
                layerstructure_NCF[-1] += 1
            else:
                layerstructure_NCF.append(1)
                current_el = el
        layerstructure_NCF[-1] += 1
        r = len(layerstructure_NCF)
    return (r, layerstructure_NCF)


## 3) Methods to randomly generate Boolean functions (uniform distribution) and Boolean networks
def random_function(n, probability_one=0.5):
    """
    Generate a random Boolean function in n variables.

    The Boolean function is represented as a truth table (an array of length 2^n) in which each entry is 0 or 1.
    Each entry is set to 1 with probability `probability_one`.

    Parameters:
        n (int): Number of variables.
        probability_one (float, optional): Probability that a given entry is 1 (default is 0.5).

    Returns:
        np.array: Boolean function as an array of length 2^n.
    """
    return np.array(np.random.random(2**n) < probability_one, dtype=int)


def random_linear_function(n):
    """
    Generate a random linear Boolean function in n variables.

    A random linear Boolean function is constructed by randomly choosing whether to include each variable or its negation in a linear sum.
    The resulting expression is then reduced modulo 2.

    Parameters:
        n (int): Number of variables.

    Returns:
        np.array: Boolean function as an array of length 2^n (truth table).
    """
    expr = '(%s) %% 2 == 1' % (' + '.join(['x%i' % i if random.random() > 0.5 else '(1 + x%i)' % i for i in range(n)]))
    return f_from_expression(expr)[0]


def random_non_degenerated_function(n, probability_one=0.5):
    """
    Generate a random non-degenerated Boolean function in n variables.

    A non-degenerated Boolean function is one in which every variable is essential (i.e. the output depends on every input).
    The function is repeatedly generated with the specified bias until a non-degenerated function is found.

    Parameters:
        n (int): Number of variables.
        probability_one (float, optional): Bias of the Boolean function (probability of a 1; default is 0.5).

    Returns:
        np.array: Boolean function as an array of length 2^n.
    
    References:
        Kadelka, C., Kuipers, J., & Laubenbacher, R. (2017). The influence of canalization on the robustness 
        of Boolean networks. Physica D: Nonlinear Phenomena, 353, 39-47.
    """
    while True:  # works because most functions are non-degenerated
        f = np.array(np.random.random(2**n) < probability_one, dtype=int)
        if not is_degenerated(f):
            return f


def random_degenerated_function(n, probability_one=0.5):
    """
    Generate a random degenerated Boolean function in n variables.

    A degenerated Boolean function is one in which at least one variable is non‐essential (its value never affects the output).
    The function is generated repeatedly until a degenerated function is found.

    Parameters:
        n (int): Number of variables.
        probability_one (float, optional): Bias of the Boolean function (default is 0.5).

    Returns:
        np.array: Boolean function as an array of length 2^n.
    
    References:
        Kadelka, C., Kuipers, J., & Laubenbacher, R. (2017). The influence of canalization on the robustness 
        of Boolean networks. Physica D: Nonlinear Phenomena, 353, 39-47.
    """
    while True:  # works not well because most functions are non-degenerated
        f = np.array(np.random.random(2**n) < probability_one, dtype=int)
        if is_degenerated(f):
            return f


def random_non_canalizing_function(n, probability_one=0.5):
    """
    Generate a random non-canalizing Boolean function in n (>1) variables.

    A Boolean function is canalizing if there exists at least one variable whose fixed value forces the output.
    This function returns one that is not canalizing.

    Parameters:
        n (int): Number of variables (n > 1).
        probability_one (float, optional): Bias of the Boolean function (default is 0.5).

    Returns:
        np.array: Boolean function as an array of length 2^n.
    
    References:
        Kadelka, C., Kuipers, J., & Laubenbacher, R. (2017). The influence of canalization on the robustness 
        of Boolean networks. Physica D: Nonlinear Phenomena, 353, 39-47.
    """
    assert n > 1
    while True:  # works because most functions are non-canalizing
        f = np.array(np.random.random(2**n) < probability_one, dtype=int)
        if not is_canalizing(f, n):
            return f


def random_non_canalizing_non_degenerated_function(n, probability_one=0.5):
    """
    Generate a random Boolean function in n (>1) variables that is both non-canalizing and non-degenerated.

    Such a function has every variable essential and is not canalizing.

    Parameters:
        n (int): Number of variables (n > 1).
        probability_one (float, optional): Bias of the Boolean function (default is 0.5).

    Returns:
        np.array: Boolean function as an array of length 2^n.
    
    References:
        Kadelka, C., Kuipers, J., & Laubenbacher, R. (2017). The influence of canalization on the robustness 
        of Boolean networks. Physica D: Nonlinear Phenomena, 353, 39-47.
    """
    assert n > 1
    while True:  # works because most functions are non-canalizing and non-degenerated
        f = np.array(np.random.random(2**n) < probability_one, dtype=int)
        if not is_canalizing(f, n) and not is_degenerated(f):
            return f


def random_k_canalizing_return_polynomial(n, k, EXACT_DEPTH_K=False, left_side_of_truth_table=[]):
    """
    Generate a random k-canalizing Boolean function in n variables and return its polynomial representation.

    A Boolean function is k-canalizing if it has at least k conditionally canalizing variables.
    The output is given as a tuple containing a list of monomials (each corresponding to a canalizing layer),
    a string representation of the core polynomial, and the primary canalized output.

    Parameters:
        n (int): Total number of variables.
        k (int): Number of canalizing variables.
        EXACT_DEPTH_K (bool, optional): If True, the generated function has exactly k canalizing variables;
                                        otherwise, the canalizing depth may exceed k (default is False).
        left_side_of_truth_table (optional): Not used in this function but provided for consistency.

    Returns:
        tuple: (monomials, core_polynomial_text, q) where:
            - monomials (list): List of strings representing monomials corresponding to the canalizing layers.
            - core_polynomial_text (str): String representation of the core polynomial.
            - q (int): The canalized output value of the first canalizing variable.
    """
    assert k > 0
    monomials = []
    aas = np.random.randint(2, size=k)  # canalizing input values
    bbs = np.random.randint(2, size=k)   # canalized output values
    if k == n and bbs[-1] != bbs[-2]:
        bbs[-1] = bbs[-2]
        aas[-1] = 1 - aas[-1]
    can_vars = np.random.choice(n, n, replace=False)
    for ii, (a, b, x_i) in enumerate(zip(aas, bbs, can_vars[:k])):
        if ii == 0:
            current_monomial = ['(x' + str(x_i) + ' + ' + str(a) + ')']
        else:
            if b == bbs[ii - 1]:
                current_monomial.append('(x' + str(x_i) + ' + ' + str(a) + ')')
            else:
                monomials.append('*'.join(current_monomial))
                current_monomial = ['(x' + str(x_i) + ' + ' + str(a) + ')']
    monomials.append('*'.join(current_monomial))
    q = bbs[0]
    if k == n:
        core_polynomial_text = '1'
        return monomials, core_polynomial_text, q
    else:
        if EXACT_DEPTH_K:
            core_polynomial = random_non_canalizing_non_degenerated_function(n - k)
        else:
            core_polynomial = random_non_degenerated_function(n - k)
        core_polynomial_text = '(' + bool_to_poly(core_polynomial, indices=can_vars[k:]) + ')'
        return monomials, core_polynomial_text, q


def random_k_canalizing(n, k, EXACT_DEPTH_K=False, left_side_of_truth_table=[], activator_or_inhibitor=[]):
    """
    Generate a random k-canalizing Boolean function in n variables.

    A Boolean function is k-canalizing if it has at least k conditionally canalizing variables.
    If EXACT_DEPTH_K is True, the function will have exactly k canalizing variables; otherwise, its canalizing depth may exceed k.

    Parameters:
        n (int): Total number of variables.
        k (int): Number of canalizing variables.
        EXACT_DEPTH_K (bool, optional): If True, enforce that the canalizing depth is exactly k (default is False).
        left_side_of_truth_table (optional): Precomputed left-hand side of the truth table for speed-up.
        activator_or_inhibitor (optional): Placeholder for future use; currently not utilized.

    Returns:
        np.array: Boolean function as an array of length 2^n (truth table).
    
    References:
        [1] He, Q., & Macauley, M. (2016). Stratification and enumeration of Boolean functions by canalizing depth. 
            Physica D: Nonlinear Phenomena, 314, 1-8.
        [2] Dimitrova, E., Stigler, B., Kadelka, C., & Murrugarra, D. (2022). Revealing the canalizing structure of Boolean functions: 
            Algorithms and applications. Automatica, 146, 110630.
    """
    try:
        assert (n - k != 1 or EXACT_DEPTH_K == False)
    except AssertionError:
        print('There are no functions of exact canalizing depth n-1.\nEither set EXACT_DEPTH_K=False or ensure k != n-1')
        return
    try:
        assert 0 <= k and k <= n
    except AssertionError:
        print('Error:\nEnsure 0 <= k <= n.')
        return
    if left_side_of_truth_table == []:  # to reduce run time, this should be computed once and then passed as argument
        left_side_of_truth_table = list(itertools.product([0, 1], repeat=n))
    num_values = 2**n
    aas = np.random.randint(2, size=k)  # canalizing inputs
    bbs = np.random.randint(2, size=k)  # canalized outputs

    # The activator_or_inhibitor parameter is currently not used.
    can_vars = np.random.choice(n, k, replace=False)
    f = np.zeros(num_values, dtype=int)
    if k < n:
        if EXACT_DEPTH_K:
            core_polynomial = random_non_canalizing_non_degenerated_function(n - k)
        else:
            core_polynomial = random_non_degenerated_function(n - k)
    else:
        core_polynomial = [1 - bbs[-1]]
    counter_non_canalized_positions = 0
    for i in range(num_values):
        for j in range(k):
            if left_side_of_truth_table[i][can_vars[j]] == aas[j]:
                f[i] = bbs[j]
                break
        else:
            f[i] = core_polynomial[counter_non_canalized_positions]
            counter_non_canalized_positions += 1
    return f


def random_k_canalizing_with_specific_layerstructure(n, layerstructure, EXACT_DEPTH_K=False, left_side_of_truth_table=[]):
    """
    Generate a random Boolean function in n variables with a specified canalizing layer structure.

    The layer structure is given as a list [k_1, ..., k_r], where each k_i indicates the number of canalizing variables 
    in that layer. If the function is fully canalizing (i.e. sum(layerstructure) == n and n > 1), the last layer must have at least 2 variables.

    Parameters:
        n (int): Total number of variables.
        layerstructure (list): List [k_1, ..., k_r] describing the canalizing layer structure.
                               Each k_i ≥ 1, and if sum(layerstructure) == n and n > 1, then layerstructure[-1] ≥ 2.
        EXACT_DEPTH_K (bool, optional): If True, enforce that the canalizing depth is exactly sum(layerstructure) (default is False).
        left_side_of_truth_table (optional): Precomputed left-hand side of the truth table for speed-up.

    Returns:
        np.array: Boolean function as an array of length 2^n (truth table).
    
    References:
        [1] He, Q., & Macauley, M. (2016). Stratification and enumeration of Boolean functions by canalizing depth.
            Physica D: Nonlinear Phenomena, 314, 1-8.
        [2] Kadelka, C., Kuipers, J., & Laubenbacher, R. (2017). The influence of canalization on the robustness 
            of Boolean networks. Physica D: Nonlinear Phenomena, 353, 39-47.
    """
    k = sum(layerstructure)  # canalizing depth
    if k == 0:
        layerstructure = [0]
    try:
        assert (n - k != 1 or EXACT_DEPTH_K == False)
    except AssertionError:
        print('Error:\nThere are no functions of exact canalizing depth n-1.\nEither set EXACT_DEPTH_K=False or ensure k=sum(layerstructure)!=n.')
        return
    try:
        assert 0 <= k and k <= n
    except AssertionError:
        print('Error:\nEnsure 0 <= k = sum(layerstructure) <= n.')
        return
    try:
        assert k < n or layerstructure[-1] > 1 or n == 1
    except AssertionError:
        print('Error:\nThe last layer of an n-canalizing function (NCF) has to have size >= 2 for n > 1.\nIf k=sum(layerstructure)=n, ensure that layerstructure[-1]>=2.')
        return
    try:
        assert min(layerstructure) >= 1
    except AssertionError:
        print('Error:\nEach layer must have at least one variable (each element of layerstructure must be >= 1).')
        return
    if left_side_of_truth_table == []:  # to decrease run time, this should be computed once and then passed as argument
        left_side_of_truth_table = list(itertools.product([0, 1], repeat=n))
    num_values = 2**n
    aas = np.random.randint(2, size=k)  # canalizing inputs
    b0 = np.random.randint(2)
    bbs = [b0] * layerstructure[0]  # canalized outputs for first layer
    for i in range(1, len(layerstructure)):
        if i % 2 == 0:
            bbs.extend([b0] * layerstructure[i])
        else:
            bbs.extend([1 - b0] * layerstructure[i])
    can_vars = np.random.choice(n, k, replace=False)
    f = np.zeros(num_values, dtype=int)
    if k < n:
        if EXACT_DEPTH_K:
            core_polynomial = random_non_canalizing_non_degenerated_function(n - k)
        else:
            core_polynomial = random_non_degenerated_function(n - k)
    else:
        core_polynomial = [1 - bbs[-1]]
    counter_non_canalized_positions = 0
    for i in range(num_values):
        for j in range(k):
            if left_side_of_truth_table[i][can_vars[j]] == aas[j]:
                f[i] = bbs[j]
                break
        else:
            f[i] = core_polynomial[counter_non_canalized_positions]
            counter_non_canalized_positions += 1
    return f


def random_adj_matrix(N, ns, NO_SELF_REGULATION=True, STRONGLY_CONNECTED=False):
    """
    Generate a random adjacency matrix for a network of N nodes.

    Each node i is assigned ns[i] outgoing edges (regulators) chosen at random.
    Optionally, self-regulation (an edge from a node to itself) can be disallowed,
    and the generated network can be forced to be strongly connected.

    Parameters:
        N (int): Number of nodes.
        ns (list or array-like): List of length N specifying the number of outgoing edges for each node.
        NO_SELF_REGULATION (bool, optional): If True, self-regulation is disallowed (default is True).
        STRONGLY_CONNECTED (bool, optional): If True, the generated network is forced to be strongly connected (default is False).

    Returns:
        tuple: (matrix, indices) where:
            - matrix (np.array): An N x N adjacency matrix with entries 0 or 1.
            - indices (list): A list of length N, where each element is an array of selected target indices for the corresponding node.
    """
    matrix = np.zeros((N, N), dtype=int)
    indices = []
    for i in range(N):
        if NO_SELF_REGULATION:
            indexes = np.random.choice(np.append(np.arange(i), np.arange(i+1, N)), ns[i], replace=False)
        else:
            indexes = np.random.choice(np.arange(N), ns[i], replace=False)
        indexes = np.sort(indexes)
        indices.append(indexes)
        for index in indexes:
            matrix[i][index] = 1
    if STRONGLY_CONNECTED:
        G = nx.from_numpy_array(matrix, parallel_edges=False, create_using=nx.MultiDiGraph())
        if not nx.is_strongly_connected(G):
            return random_adj_matrix(N, ns, NO_SELF_REGULATION, STRONGLY_CONNECTED)
    return (matrix, indices)


def get_strongly_connected_components(I):
    """
    Compute the strongly connected components of a wiring diagram.

    The wiring diagram is provided as a list of lists I, where I[i] contains the indices of regulators for node i.
    The function constructs a directed graph from these edges and returns its strongly connected components.

    Parameters:
        I (list): A list of lists, where each inner list contains the regulators (source nodes) for the corresponding target node.

    Returns:
        list: A list of sets, each representing a strongly connected component.
    """
    edges_wiring_diagram = []
    for target, regulators in enumerate(I):
        for regulator in regulators:
            edges_wiring_diagram.append((regulator, target))
    subG = nx.from_edgelist(edges_wiring_diagram, create_using=nx.MultiDiGraph())
    return [scc for scc in nx.strongly_connected_components(subG)]


def random_edge_list(N, ns, NO_SELF_REGULATION, AT_LEAST_ONE_REGULATOR_PER_GENE=False):
    """
    Generate a random edge list for a network of N nodes with optional constraints.

    Each node i receives ns[i] incoming edges chosen at random.
    Optionally, the function ensures that every node has at least one regulator.

    Parameters:
        N (int): Number of nodes.
        ns (list or array-like): List of length N specifying the number of regulators for each node.
        NO_SELF_REGULATION (bool): If True, disallow self-regulation.
        AT_LEAST_ONE_REGULATOR_PER_GENE (bool, optional): If True, ensure that each node has at least one outgoing edge (default is False).

    Returns:
        list: A list of tuples (source, target) representing the edges.
    """
    if AT_LEAST_ONE_REGULATOR_PER_GENE == False:
        edge_list = []
        for i in range(N):
            if NO_SELF_REGULATION:
                indices = np.random.choice(np.append(np.arange(i), np.arange(i+1, N)), ns[i], replace=False)
            else:
                indices = np.random.choice(np.arange(N), ns[i], replace=False)
            edge_list.extend(list(zip(indices, i * np.ones(ns[i], dtype=int))))
    else:
        edge_list = []
        outdegree = np.zeros(N, dtype=int)
        sum_ns = sum(ns)  # total number of regulations
        for i in range(N):
            if NO_SELF_REGULATION:
                indices = np.random.choice(np.append(np.arange(i), np.arange(i+1, N)), ns[i], replace=False)
            else:
                indices = np.random.choice(np.arange(N), ns[i], replace=False)
            outdegree[indices] += 1
            edge_list.extend(list(zip(indices, i * np.ones(ns[i], dtype=int))))
        while min(outdegree) == 0:
            index_sink = np.where(outdegree == 0)[0][0]
            index_edge = int(random.random() * sum_ns)
            if NO_SELF_REGULATION:
                while edge_list[index_edge][1] == index_sink:
                    index_edge = int(random.random() * sum_ns)
            outdegree[index_sink] += 1
            outdegree[edge_list[index_edge][0]] -= 1
            edge_list[index_edge] = (index_sink, edge_list[index_edge][1])
    return edge_list

def get_essential_network(F, I):
    """
    Determine the essential components of a Boolean network.

    For each gene in a Boolean network, represented by its Boolean function and its regulators,
    this function extracts the “essential” part of the function by removing non-essential regulators.
    The resulting network contains, for each gene, a reduced truth table (with only the essential inputs)
    and a corresponding list of essential regulators.

    Parameters:
        F (list): A list of N Boolean functions (truth tables). For gene i, the Boolean function is given as a list
                  of length 2^(n_i), where n_i is the number of regulators for that gene.
        I (list): A list of N lists. For gene i, I[i] is a list of regulator indices (typically 0, 1, ..., n_i-1)
                  corresponding to the wiring diagram of the Boolean network.

    Returns:
        tuple: (F_essential, I_essential) where:
            - F_essential is a list of N Boolean functions (truth tables) of length 2^(m_i), with m_i ≤ n_i,
              representing the functions restricted to the essential regulators.
            - I_essential is a list of N lists containing the indices of the essential regulators for each gene.
    """
    import itertools
    F_essential = []
    I_essential = []
    for f, regulators in zip(F, I):
        if len(f) == 0:  # happens if the actual degree of f was too large for it to be loaded
            F_essential.append(f)
            I_essential.append(regulators)
            continue
        elif sum(f) == 0:
            F_essential.append(np.array([0]))
            I_essential.append(np.array([], dtype=int))
            continue
        elif sum(f) == len(f):
            F_essential.append(np.array([1]))
            I_essential.append(np.array([], dtype=int))
            continue
        essential_variables = np.array(get_essential_variables(f))
        n = len(regulators)
        non_essential_variables = np.array(list(set(list(range(n))) - set(essential_variables)))
        if len(non_essential_variables) == 0:
            F_essential.append(f)
            I_essential.append(regulators)
        else:
            left_side_of_truth_table = np.array(list(itertools.product([0, 1], repeat=n)))
            F_essential.append(np.array(f)[np.sum(left_side_of_truth_table[:, non_essential_variables], 1) == 0])
            I_essential.append(np.array(regulators)[essential_variables])
    return F_essential, I_essential


def random_BN(N, n=2, k=0, STRONGLY_CONNECTED=True, indegree_distribution='constant',
              list_x=[], kis=None, EXACT_DEPTH=False, NO_SELF_REGULATION=True, LINEAR=False,
              edges_wiring_diagram=None, bias=0.5):
    """
    Generate a random Boolean network (BN).

    This function creates a Boolean network for N genes by first generating a wiring diagram
    (a set of regulatory interactions) according to a specified in-degree distribution and then assigning
    Boolean functions to each gene. The functions can be generated with prescribed canalizing properties,
    linearity, and bias. Optionally, user-defined in-degrees, canalizing depths (k or kis), or even an edge list
    (wiring diagram) can be provided.

    Parameters:
        N (int): Number of genes (nodes) in the network.
        n (int, list, or np.array, optional): Determines the in-degree for each gene. If an integer, every gene has the same number of regulators;
                                                if a vector, each element gives the number of regulators for the corresponding gene.
                                                Default is 2.
        k (int, list, or np.array, optional): Specifies the canalizing depth for each gene.
                                              If an integer, the same depth is used for all genes; if a vector, each gene gets its own depth.
                                              Default is 0.
        STRONGLY_CONNECTED (bool, optional): If True, ensures that the generated network is strongly connected. Default is True.
        indegree_distribution (str, optional): In-degree distribution to use. Options include 'constant' (or 'dirac'/'delta'),
                                               'uniform', or 'poisson'. Default is 'constant'.
        list_x (list, optional): Precomputed truth tables (lists of tuples) for different in-degrees, used for speed-up. Default is an empty list.
        kis (optional): Specifies the canalizing layer structure for the Boolean functions. If provided, the parameter k is ignored.
        EXACT_DEPTH (bool, optional): If True, Boolean functions are generated with exactly the specified canalizing depth;
                                      if False, the functions have at least that depth. Default is False.
        NO_SELF_REGULATION (bool, optional): If True, self-regulation (self-loops) is disallowed. Default is True.
        LINEAR (bool, optional): If True, Boolean functions are generated to be linear. Default is False.
        edges_wiring_diagram (optional): User-defined edge list for the wiring diagram. If provided, this is used instead of generating one.
        bias (float, optional): Bias for generating Boolean functions (probability of output 1). Default is 0.5.

    Returns:
        tuple: (F, I, ns) where:
            - F is a list of N Boolean functions (truth tables) for the genes.
            - I is the wiring diagram, given as a list of N lists; I[i] contains the indices of regulators for gene i.
            - ns is an array (or list) of length N representing the in-degree (number of regulators) for each gene.
    """
    # Process the input for in-degree distribution
    if indegree_distribution in ['constant', 'dirac', 'delta']:
        if type(n) in [list, np.array]:
            try:
                assert type(n) in [list, np.array]
                assert np.all([type(el) in [int, np.int64] for el in n])
                assert len(n) == N
                assert min(n) >= 1
                assert max(n) <= N
                ns = np.array(n[:])
            except AssertionError:
                print('Error: A vector n was submitted.\nTo use a user-defined in-degree vector, ensure that n is an N-dimensional vector where each element is an integer between 1 and N.')
                return
        else:
            try:
                assert type(n) in [int, np.int64]
                assert n >= 1
                assert n <= N
                ns = np.ones(N, dtype=int) * n
            except AssertionError:
                print('Error: n must be a single integer (or N-dimensional vector of integers) between 1 and N when using a constant degree distribution.')
                return
    elif indegree_distribution == 'uniform':
        if type(n) in [list, np.array]:
            try:
                assert type(n) in [list, np.array]
                assert np.all([type(el) in [int, np.int64] for el in n])
                assert len(n) == N
                assert min(n) >= 1
                assert max(n) <= N
                ns = np.array(n[:])
            except AssertionError:
                print('Error: A vector n was submitted.\nEnsure that n is an N-dimensional vector where each element is an integer between 1 and N representing the upper bound of a uniform degree distribution (lower bound == 1).')
                return
        else:
            try:
                assert type(n) in [int, np.int64]
                assert n >= 1
                assert n <= N
                ns = np.ones(N, dtype=int) * n
            except AssertionError:
                print('Error: n must be a single integer (or N-dimensional vector of integers) between 1 and N representing the upper bound of a uniform degree distribution (lower bound == 1).')
                return
    elif indegree_distribution == 'poisson':
        if type(n) in [list, np.array]:
            try:
                assert type(n) in [list, np.array]
                assert np.all([type(el) in [int, np.int64, float, np.float64] for el in n])
                assert len(n) == N
                assert min(n) >= 1
                assert max(n) <= N
                ns = np.array(n[:])
            except AssertionError:
                print('Error: A vector n was submitted.\nEnsure that n is an N-dimensional vector where each element is > 0 and represents the Poisson parameter.')
                return
        else:
            try:
                assert type(n) in [int, np.int64, float, np.float64]
                assert n >= 1
                assert n <= N
                ns = np.ones(N, dtype=int) * n
            except AssertionError:
                print('Error: n must be a single number (or N-dimensional vector) > 0 representing the Poisson parameter.')
                return
    else:
        print('None of the predefined in-degree distributions were chosen.\nTo use a user-defined in-degree vector, use the input n to submit an N-dimensional vector where each element of n must be between 1 and N.')
        return

    if kis is None:
        if type(k) in [int, np.int64]:
            try:
                assert k >= 0
                assert k <= N
                max_k = k
            except AssertionError:
                print('Error: k must be an integer between 0 and N.')
                return
        elif type(k) in [list, np.array]:
            try:
                assert len(k) == N
                assert np.all([type(el) in [int, np.int64] for el in k])
                max_k = max(k)
                assert min(k) >= 0
                assert max_k <= N
            except AssertionError:
                print('Error: A vector k was submitted.\nTo use a user-defined vector k, ensure that k is an N-dimensional vector where each element is an integer between 0 and N.')
                return
        else:
            print('Error: Wrong input format for k.\nk must be a single integer (or N-dimensional vector of integers) between 0 and N.')
            return
    else:  # kis provided
        if np.all([type(el) in [int, np.int64] for el in kis]):
            try:
                assert np.all([type(el) in [int, np.int64] for el in kis])
                assert np.all([el >= 1 for el in kis])
                max_k = sum(kis)
                assert max_k <= N
            except AssertionError:
                print('Error: The layer structure kis must be a vector of positive integers with 0 <= k = sum(kis) <= N.')
                return
        elif np.all([type(el) in [list, np.array] for el in kis]):
            try:
                assert len(kis) == N
                assert type(kis[0][0]) in [int, np.int64]
                max_k = max([sum(el) for el in kis])
                assert min([min(el) for el in kis]) >= 0
                assert max_k <= N
            except AssertionError:
                print('Error: A vector of kis was submitted.\nEnsure that kis is an N-dimensional vector where each element represents a layer structure and is a vector of positive integers with 1 <= k = sum(kis[i]) <= N.')
                return
        else:
            print('Error: Wrong input format for kis.\nkis must be a single vector (or N-dimensional vector of layer structures) where the sum of each element is between 0 and N.')
            return

    if edges_wiring_diagram is None:
        while True:  # Keep generating until we have a strongly connected graph
            if indegree_distribution == 'uniform':
                ns = 1 + np.random.randint(n - 1, size=N)
            elif indegree_distribution == 'poisson':
                ns = np.random.poisson(lam=n, size=N)
                ns[ns == 0] = 1
                ns[ns > N - int(NO_SELF_REGULATION)] = N - int(NO_SELF_REGULATION)
            edges_wiring_diagram = random_edge_list(N, ns, NO_SELF_REGULATION)
            if STRONGLY_CONNECTED:
                G = nx.from_edgelist(edges_wiring_diagram, create_using=nx.MultiDiGraph())
                if not nx.is_strongly_connected(G):
                    continue
            break
    else:
        try:
            assert len(set(np.array(edges_wiring_diagram).flatten())) == N
        except AssertionError:
            print("Number of nodes provided in edges_wiring_diagram != N")
            return
        ns = np.zeros(N, dtype=int)
        for target in np.array(edges_wiring_diagram)[:, 1]:
            ns[target] += 1

    max_n = max(ns)
    if max_k > 0 and (list_x == [] or len(list_x) < max_n):
        list_x = [[[0], [1]]]
        list_x.extend([list(itertools.product([0, 1], repeat=nn)) for nn in range(2, max_n + 1)])

    F = []
    for i in range(N):
        if LINEAR:
            F.append(random_linear_function(ns[i]))
        if k > 0 and kis is None:
            if type(k) in [int, np.int64]:
                F.append(random_k_canalizing(ns[i], min(k, ns[i]), EXACT_DEPTH_K=EXACT_DEPTH, left_side_of_truth_table=list_x[ns[i]-1]))
            else:
                F.append(random_k_canalizing(ns[i], min(k[i], ns[i]), EXACT_DEPTH_K=EXACT_DEPTH, left_side_of_truth_table=list_x[ns[i]-1]))
        elif kis is not None:
            if np.all([type(el) in [int, np.int64] for el in kis]):
                F.append(random_k_canalizing_with_specific_layerstructure(ns[i], kis, EXACT_DEPTH_K=EXACT_DEPTH, left_side_of_truth_table=list_x[ns[i]-1]))
            else:
                F.append(random_k_canalizing_with_specific_layerstructure(ns[i], kis[i], EXACT_DEPTH_K=EXACT_DEPTH, left_side_of_truth_table=list_x[ns[i]-1]))
        else:
            if EXACT_DEPTH is True:
                F.append(random_non_canalizing_non_degenerated_function(ns[i], bias))
            else:
                F.append(random_non_degenerated_function(ns[i], bias))

    I = [[] for _ in range(N)]
    for edge in edges_wiring_diagram:
        I[edge[1]].append(edge[0])
    for i in range(N):
        I[i] = np.sort(I[i])
    return F, I, ns


def get_perturbed_network(F, I, ns, control_target, control_source, type_of_control=0, left_side_of_truth_table=[]):
    """
    Generate a perturbed Boolean network by removing the influence of a specified regulator.

    The function modifies the Boolean function for a target gene by restricting it to those entries in its truth table
    where the input from a given regulator equals the specified type_of_control. The regulator is then removed from
    the wiring diagram for that gene.

    Parameters:
        F (list): List of Boolean functions (truth tables) for each gene.
        I (list): Wiring diagram for the network; each entry I[i] is a list of regulator indices for gene i.
        ns (list or np.array): List of in-degrees (number of regulators) for each gene.
        control_target (int): Index of the target gene to be perturbed.
        control_source (int): Index of the regulator whose influence is to be removed.
        type_of_control (int, optional): The regulator value (0 or 1) for which the perturbation is applied. Default is 0.
        left_side_of_truth_table (optional): Precomputed truth table (array of tuples) for the target gene with ns[control_target] inputs.
                                              If not provided, it is computed.

    Returns:
        tuple: (F_new, I_new, ns) where:
            - F_new is the updated list of Boolean functions after perturbation.
            - I_new is the updated wiring diagram after removing the control regulator from the target gene.
            - ns is the updated list of in-degrees for each gene.
    """
    F_new = [f for f in F]
    I_new = [i for i in I]

    if left_side_of_truth_table == []:
        left_side_of_truth_table = np.array(list(itertools.product([0, 1], repeat=ns[control_target])))

    try:
        index = list(I[control_target]).index(control_source)
        F_new[control_target] = F_new[control_target][left_side_of_truth_table[:, index] == type_of_control]
        dummy = list(I_new[control_target])
        dummy.remove(control_source)
        I_new[control_target] = np.array(dummy)
    except ValueError:
        print('source not regulating target')

    ns = list(map(len, I_new))
    return F_new, I_new, ns

def get_constant_nodes(I, degree, N):
    """
    Identify constant nodes in a Boolean network.

    A node is considered constant if it has exactly one regulator and that regulator is the node itself.

    Parameters:
        I (list): A list where I[i] is a list of regulator indices for node i.
        degree (list or array): A list where degree[i] is the number of regulators for node i.
        N (int): Total number of nodes in the network.

    Returns:
        np.array: Array of node indices that are constant.
    """
    return np.array([i for i in range(N) if degree[i] == 1 and I[i][0] == i])

## 5) Analysis methods
def update_single_node(f, states_regulators):
    """
    Update the state of a single node.

    The new state is obtained by applying the Boolean function f to the states of its regulators.
    The regulator states are converted to a decimal index using bin2dec.

    Parameters:
        f (list or np.array): Boolean function (truth table) for the node.
        states_regulators (list or np.array): Binary vector representing the states of the node's regulators.

    Returns:
        int: Updated state of the node (0 or 1).
    """
    return f[bin2dec(states_regulators)]

def update_network_synchronously(F, I, X):
    """
    Perform a synchronous update of a Boolean network.

    Each node's new state is determined by applying its Boolean function to the current states of its regulators.
    The conversion from the regulator states (a binary vector) to a truth table index is done via bin2dec.

    Parameters:
        F (list): List of Boolean functions (truth tables) for each node.
        I (list): List of regulator indices for each node.
        X (list or np.array): Current state vector of the network.

    Returns:
        np.array: New state vector after the update.
    """
    if type(X)==list:
        X = np.array(X)
    Fx = np.zeros(len(F), dtype=int)
    for i in range(len(F)):
        Fx[i] = update_single_node(f = F[i], states_regulators = X[I[i]])
    return Fx


def update_network_SDDS(F, I, N, X, P):
    """
    Perform a stochastic update (SDDS) on a Boolean network.

    For each node, the next state is computed as nextstep = F[i] evaluated on the current states of its regulators.
    If nextstep > X[i], the node is activated with probability P[i,0]; if nextstep < X[i],
    the node is degraded with probability P[i,1]. Otherwise, the state remains unchanged.

    Parameters:
        F (list): List of Boolean functions (truth tables) for each node.
        I (list): List of regulator indices for each node.
        X (list or np.array): Current state vector.
        P (np.array): A len(F)×2 array of probabilities; for each node i, P[i,0] is the activation probability,
                      and P[i,1] is the degradation probability.

    Returns:
        np.array: Updated state vector after applying the stochastic update.
    """
    if type(X)==list:
        X = np.array(X)
    Fx = X.copy()
    for i in range(len(F)):
        nextstep = update_single_node(f = F[i], states_regulators = X[I[i]])
        if nextstep > X[i] and random.random() < P[i, 0]:  # activation
            Fx[i] = nextstep
        elif nextstep < X[i] and random.random() < P[i, 1]:  # degradation
            Fx[i] = nextstep
    return Fx


def update_network_many_times(F, I, initial_state, n_steps):
    """
    Update the state of a Boolean network sychronously multiple time steps.

    Starting from the initial state, the network is updated synchronously n_steps times using the update_network_synchronously function.

    Parameters:
        F (list): List of Boolean functions for each node.
        I (list): List of regulator indices for each node.
        initial_state (list or np.array): Initial state vector of the network.
        n_steps (int): Number of update iterations to perform.

    Returns:
        np.array: Final state vector after n_steps updates.
    """
    N = len(F)
    for i in range(n_steps):
        initial_state = update_network_synchronously(F, I, N, initial_state)
    return initial_state


def get_derrida_value(F, I, nsim):
    """
    Estimate the Derrida value for a Boolean network.

    The Derrida value is computed by perturbing a single node in a randomly chosen state and measuring
    the average Hamming distance between the resulting updated states of the original and perturbed networks.

    Parameters:
        F (list): List of Boolean functions (truth tables) for each node.
        I (list): List of regulator indices for each node.
        nsim (int): Number of simulations to perform.

    Returns:
        float: The average Hamming distance (Derrida value) over nsim simulations.
    """
    N = len(F)
    hamming_distances = []
    for i in range(nsim):
        X = np.random.randint(0, 2, N)
        Y = X.copy()
        index = np.random.randint(N)
        Y[index] = 1 - Y[index]
        FX = update_network_synchronously(F, I, X)
        FY = update_network_synchronously(F, I, Y)
        hamming_distances.append(sum(FX != FY))
    return np.mean(hamming_distances)



def get_steady_states_asynchronous(F, I, N, nsim=500, EXACT=False, left_side_of_truth_table=[], 
                                   initial_sample_points=[], search_depth=50, SEED=-1, DEBUG=True):
    """
    Compute the steady states of a Boolean network under asynchronous updates.

    This function simulates asynchronous updates of a Boolean network (with N nodes)
    for a given number of initial conditions (nsim). For each initial state, the network
    is updated asynchronously until a steady state (or attractor) is reached or until a maximum
    search depth is exceeded. The simulation can be performed either approximately (by sampling nsim
    random initial conditions) or exactly (by iterating over the entire state space when EXACT=True).

    Parameters:
        F (list): List of Boolean functions (truth tables) for each node.
                 Each function is defined over 2^(# of regulators) entries.
        I (list): List of lists, where I[i] contains the indices of the regulators for node i.
        N (int): Total number of nodes in the network.
        nsim (int, optional): Number of initial conditions to simulate (default is 500).
        EXACT (bool, optional): If True, iterate over the entire state space (2^N initial conditions);
                                otherwise, use nsim random initial conditions. (Default is False.)
        left_side_of_truth_table (list, optional): Precomputed truth table (list of tuples) for N inputs.
                                                     Used only if EXACT is True.
        initial_sample_points (list, optional): List of initial states (as binary vectors) to use.
                                                  If provided and EXACT is False, these override random sampling.
        search_depth (int, optional): Maximum number of asynchronous update iterations to attempt per simulation.
        SEED (int, optional): Random seed. If SEED is -1, a random seed is generated.
        DEBUG (bool, optional): If True, print debugging information during simulation.

    Returns:
        tuple: A tuple containing:
            - steady_states (list): List of steady state values (in decimal form) found.
            - number_of_steady_states (int): Total number of unique steady states.
            - basin_sizes (list): List of counts showing how many initial conditions converged to each steady state.
            - steady_state_dict (dict): Dictionary mapping a steady state (in decimal) to its index in the steady_states list.
            - dictF (dict): Dictionary caching state transitions. Keys are tuples (xdec, i) and values are the updated state.
            - SEED (int): The random seed used for the simulation.
            - initial_sample_points (list): The list of initial sample points used (if provided) or those generated during simulation.
    """
    if EXACT and left_side_of_truth_table == []:
        left_side_of_truth_table = list(map(np.array, list(itertools.product([0, 1], repeat=N))))

    sampled_points = []
    
    assert initial_sample_points == [] or not EXACT, (
        "Warning: sample points were provided but, with option EXACT==True, the entire state space is computed "
        "(and initial sample points ignored)"
    )
    
    if SEED == -1:
        SEED = int(random.random() * 2**31)
    
    np.random.seed(SEED)
    
    dictF = dict()
    steady_states = []
    basin_sizes = []
    steady_state_dict = dict()   
    
    for iteration in range(nsim if not EXACT else 2**N):
        if EXACT:
            x = left_side_of_truth_table[iteration]
            xdec = iteration
        else:
            if initial_sample_points == []:  # generate random initial states on the fly
                x = np.random.randint(2, size=N)
                xdec = bin2dec(x)
                sampled_points.append(xdec)
            else:                
                x = initial_sample_points[iteration]
                xdec = bin2dec(x)
        
        if DEBUG:
            print(iteration, -1, -1, False, xdec, x)
        for jj in range(search_depth):  # update until a steady state is reached or search_depth is exceeded
            FOUND_NEW_STATE = False
            try:
                # Check if this state is already recognized as a steady state.
                index_ss = steady_state_dict[xdec]
            except KeyError:
                # Asynchronously update the state until a new state is found.
                update_order_to_try = np.random.permutation(N)
                for i in update_order_to_try:
                    try:
                        fxdec = dictF[(xdec, i)]
                        if fxdec != xdec:
                            FOUND_NEW_STATE = True
                            x[i] = 1 - x[i]
                    except KeyError:
                        fx_i = update_single_node(F[i], x[I[i]])
                        if fx_i > x[i]:
                            fxdec = xdec + 2**(N - 1 - i)
                            x[i] = 1
                            FOUND_NEW_STATE = True
                        elif fx_i < x[i]:
                            fxdec = xdec - 2**(N - 1 - i)
                            x[i] = 0
                            FOUND_NEW_STATE = True
                        else:
                            fxdec = xdec
                        dictF.update({(xdec, i): fxdec})
                    if FOUND_NEW_STATE:
                        xdec = fxdec
                        break
                if DEBUG:
                    print(iteration, jj, i, FOUND_NEW_STATE, xdec, x)
            if FOUND_NEW_STATE == False:  # steady state reached
                try:
                    index_ss = steady_state_dict[xdec]
                    basin_sizes[index_ss] += 1
                    break
                except KeyError:
                    steady_state_dict.update({xdec: len(steady_states)})
                    steady_states.append(xdec)
                    basin_sizes.append(1)
                    break
        if DEBUG:
            print()
    if sum(basin_sizes) < (nsim if not EXACT else 2**N):
        print('Warning: only %i of the %i tested initial conditions eventually reached a steady state. Try increasing the search depth. '
              'It may however also be the case that your asynchronous state space contains a limit cycle.' %
              (sum(basin_sizes), nsim if not EXACT else 2**N))
    return (steady_states, len(steady_states), basin_sizes, steady_state_dict, dictF, SEED,
            initial_sample_points if initial_sample_points != [] else sampled_points)


def get_steady_states_asynchronous_given_one_initial_condition(F, I, nsim=500, stochastic_weights=[], initial_condition=0, search_depth=50, SEED=-1, DEBUG=False):
    """
    Determine the steady states reachable from one initial condition using weighted asynchronous updates.

    This function is similar to steady_states_asynchronous_given_one_IC but allows the update order
    to be influenced by provided stochastic weights (one per node). A weight vector (of length N) may be provided,
    and if given, it is normalized and used to bias the random permutation of node update order.
    
    Parameters:
        F (list): List of Boolean functions (truth tables) for each node.
        I (list): List of regulator indices for each node.
        nsim (int, optional): Number of simulation runs (default is 500).
        stochastic_weights (list, optional): List of stochastic weights (one per node) used to bias update order.
                                               If empty, uniform random order is used.
        initial_condition (int or list/np.array, optional): The initial state for all simulations. If an integer, 
                                                           it is converted to a binary vector. Default is 0.
        search_depth (int, optional): Maximum number of asynchronous update iterations per simulation (default is 50).
        SEED (int, optional): Random seed. If -1, a random seed is generated (default is -1).
        DEBUG (bool, optional): If True, print debugging information (default is False).

    Returns:
        tuple: A tuple containing:
            - steady_states (list): List of steady state values (in decimal form) reached.
            - number_of_steady_states (int): Total number of unique steady states.
            - basin_sizes (list): List of counts of how many simulations reached each steady state.
            - transient_times (list): List of lists with transient times (number of updates) for each steady state.
            - steady_state_dict (dict): Dictionary mapping a steady state (in decimal) to its index.
            - dictF (dict): Dictionary caching computed state transitions.
            - SEED (int): The random seed used.
            - queues (list): List of state update queues (the sequence of states encountered) for each simulation.
    """
    if SEED == -1:
        SEED = int(random.random() * 2**31)
    np.random.seed(SEED)    
    
    N = len(F)
    
    if type(initial_condition) == int:
        initial_condition = np.array(dec2bin(initial_condition, N))
        initial_condition_bin = bin2dec(initial_condition)
    else:
        initial_condition = np.array(initial_condition, dtype=int)
        initial_condition_bin = bin2dec(initial_condition)
    
    assert stochastic_weights == [] or len(stochastic_weights) == N, "one stochastic weight per node is required"    
    if stochastic_weights != []:
        stochastic_weights = np.array(stochastic_weights) / sum(stochastic_weights)
    
    dictF = dict()
    steady_states = []
    basin_sizes = []
    transient_times = []
    steady_state_dict = dict()   
    queues = []
    for iteration in range(nsim):
        x = initial_condition.copy()
        xdec = initial_condition_bin
        queue = [xdec]
        for jj in range(search_depth):  # update until a steady state is reached or search_depth is exceeded
            FOUND_NEW_STATE = False
            try:
                index_ss = steady_state_dict[xdec]
            except KeyError:
                if stochastic_weights != []:
                    update_order_to_try = np.random.choice(N, size=N, replace=False, p=stochastic_weights)
                else:
                    update_order_to_try = np.random.permutation(N)
                for i in update_order_to_try:
                    try:
                        fxdec = dictF[(xdec, i)]
                        if fxdec != xdec:
                            FOUND_NEW_STATE = True
                            x[i] = 1 - x[i]
                    except KeyError:
                        fx_i = update_single_node(F[i], x[I[i]])
                        if fx_i > x[i]:
                            fxdec = xdec + 2**(N - 1 - i)
                            x[i] = 1
                            FOUND_NEW_STATE = True
                        elif fx_i < x[i]:
                            fxdec = xdec - 2**(N - 1 - i)
                            x[i] = 0
                            FOUND_NEW_STATE = True
                        else:
                            fxdec = xdec
                        dictF.update({(xdec, i): fxdec})
                    if FOUND_NEW_STATE:
                        xdec = fxdec
                        queue.append(xdec)
                        break
                if DEBUG:
                    print(iteration, jj, i, FOUND_NEW_STATE, xdec, x)
            if not FOUND_NEW_STATE:  # steady state reached
                queues.append(queue[:])
                try:
                    index_ss = steady_state_dict[xdec]
                    basin_sizes[index_ss] += 1
                    transient_times[index_ss].append(jj)
                    break
                except KeyError:
                    steady_state_dict.update({xdec: len(steady_states)})
                    steady_states.append(xdec)
                    basin_sizes.append(1)
                    transient_times.append([jj])
                    break
        if FOUND_NEW_STATE:
            print(jj)
            break
        if DEBUG:
            print()
    if sum(basin_sizes) < nsim:
        print('Warning: only %i of the %i tested initial conditions eventually reached a steady state. '
              'Try increasing the search depth. It may also be that your asynchronous state space contains a limit cycle.' % (sum(basin_sizes), nsim))
    return (steady_states, len(steady_states), basin_sizes, transient_times, steady_state_dict, dictF, SEED, queues)


def get_attractors_synchronous(F, I, nsim=500, initial_sample_points=[], n_steps_timeout=1000000000000,
                               INITIAL_SAMPLE_POINTS_AS_BINARY_VECTORS=True):
    """
    Compute the number of attractors in a Boolean network using an alternative (v2) approach.

    This version is optimized for networks with longer average path lengths. For each of nb initial conditions,
    the network is updated synchronously until an attractor is reached or until n_steps_timeout is exceeded.
    The function returns the attractors found, their basin sizes, a mapping of states to attractors,
    the set of initial sample points used, the explored state space, and the number of simulations that timed out.

    Parameters:
        F (list): List of Boolean functions (truth tables) for each node.
        I (list): List of lists, where I[i] contains the indices of the regulators for node i.
        nsim (int, optional): Number of initial conditions to simulate (default is 500).
        initial_sample_points (list, optional): List of initial states (in decimal) to use.
        n_steps_timeout (int, optional): Maximum number of update steps allowed per simulation (default is a very large number).
        INITIAL_SAMPLE_POINTS_AS_BINARY_VECTORS (bool, optional): If True, initial_sample_points are provided as binary vectors;
                                                                  if False, they are given as decimal numbers. Default is True.

    Returns:
        tuple: A tuple containing:
            - attractors (list): List of attractors (each as a list of states in the attractor cycle).
            - number_of_attractors (int): Total number of unique attractors found.
            - basin_sizes (list): List of counts for each attractor.
            - attr_dict (dict): Dictionary mapping states (in decimal) to the index of their attractor.
            - initial_sample_points (list): The initial sample points used (if provided, they are returned; otherwise, the generated points).
            - state_space (list): List of states (in decimal) encountered after one update from initial_sample_points.
            - n_timeout (int): Number of simulations that reached the step timeout.
    """
    dictF = dict()
    attractors = []
    basin_sizes = []
    attr_dict = dict()
    state_space = []
    
    sampled_points = []
    n_timeout = 0
    
    N = len(F)
    
    INITIAL_SAMPLE_POINTS_EMPTY = check_if_empty(initial_sample_points)
    if not INITIAL_SAMPLE_POINTS_EMPTY:
        nsim = len(initial_sample_points)
    
    for i in range(nsim):
        if INITIAL_SAMPLE_POINTS_EMPTY:
            x = np.random.randint(2, size=N)
            xdec = bin2dec(x)
            sampled_points.append(xdec)
        else:
            if INITIAL_SAMPLE_POINTS_AS_BINARY_VECTORS:
                x = initial_sample_points[i]
                xdec = bin2dec(x)
            else:
                xdec = initial_sample_points[i]
                x = np.array(dec2bin(xdec, N))
        queue = [xdec]
        count = 0
        while count < n_steps_timeout:
            try:
                fxdec = dictF[xdec]
            except KeyError:
                fx = update_network_synchronously(F, I, x)
                fxdec = bin2dec(fx)
                dictF.update({xdec: fxdec})
                x = fx
            if count == 0:
                state_space.append(fxdec)
            try:
                index_attr = attr_dict[fxdec]
                basin_sizes[index_attr] += 1
                attr_dict.update(list(zip(queue, [index_attr] * len(queue))))
                break
            except KeyError:
                try:
                    index = queue.index(fxdec)
                    attr_dict.update(list(zip(queue[index:], [len(attractors)] * (len(queue) - index))))
                    attractors.append(queue[index:])
                    basin_sizes.append(1)
                    break
                except ValueError:
                    pass
            queue.append(fxdec)
            xdec = fxdec
            count += 1
            if count == n_steps_timeout:
                n_timeout += 1            
    return (attractors, len(attractors), basin_sizes, attr_dict,
            sampled_points if INITIAL_SAMPLE_POINTS_EMPTY else initial_sample_points,
            state_space, n_timeout)


def get_attractors_synchronous_exact(F, I, left_side_of_truth_table=None):
    """
    Compute the exact number of attractors in a Boolean network using a fast, vectorized approach.

    This function computes the state of each node for all 2^N states by constructing the network's state space,
    then maps each state to its corresponding successor state via the Boolean functions F.
    Attractors and their basin sizes are then determined by iterating over the entire state space.

    Parameters:
        F (list): List of Boolean functions (truth tables) for each node.
        I (list): List of lists, where I[i] contains the indices of the regulators for node i.
        left_side_of_truth_table (np.array, optional): Precomputed array of all 2^N states (each row is a state).
                                                        If None, it is generated.

    Returns:
        tuple: A tuple containing:
            - attractors (list): List of attractors (each attractor is represented as a list of states forming the cycle).
            - number_of_attractors (int): Total number of unique attractors.
            - basin_sizes (list): List of counts for each attractor.
            - attractor_dict (dict): Dictionary mapping each state (in decimal) to its attractor index.
            - state_space (np.array): The constructed state space matrix (of shape (2^N, N)).
    """
    
    N = len(F)
    
    if left_side_of_truth_table is None:
        left_side_of_truth_table = np.array(list(map(np.array, list(itertools.product([0, 1], repeat=N)))))
    
    b_for_bin2dec = np.array([2**i for i in range(N)])[::-1]
    degrees = list(map(len, I))
    
    state_space = np.zeros((2**N, N), dtype=int)
    for i in range(N):
        for j, x in enumerate(itertools.product([0, 1], repeat=degrees[i])):
            if F[i][j]:
                # For rows in left_side_of_truth_table where the columns I[i] equal x, set state_space accordingly.
                state_space[np.all(left_side_of_truth_table[:, I[i]] == np.array(x), axis=1), i] = 1
                
    dictF = dict(zip(list(range(2**N)), np.dot(state_space, b_for_bin2dec)))
    
    attractors = []
    basin_sizes = []
    attractor_dict = dict()
    for xdec in range(2**N):
        queue = [xdec]
        while True:
            fxdec = dictF[xdec]
            try:
                index_attr = attractor_dict[fxdec]
                basin_sizes[index_attr] += 1
                attractor_dict.update(list(zip(queue, [index_attr] * len(queue))))
                break
            except KeyError:
                try:
                    index = queue.index(fxdec)
                    attractor_dict.update(list(zip(queue, [len(attractors)] * len(queue))))
                    attractors.append(queue[index:])
                    basin_sizes.append(1)
                    break
                except ValueError:
                    pass
            queue.append(fxdec)
            xdec = fxdec
    return (attractors, len(attractors), basin_sizes, attractor_dict, state_space)


def get_proportion_of_largest_basin_size(basin_sizes):
    """
    Compute the proportion of the largest basin size relative to the total basin sizes.

    This function calculates the ratio of the largest basin size to the sum of all basin sizes.
    This metric is useful for assessing the dominance of a particular attractor’s basin in a Boolean network.

    Parameters:
        basin_sizes (list or array-like): A list where each element represents the size of a basin
                                           (i.e., the number of initial conditions that converge to a specific attractor).

    Returns:
        float: The proportion of the largest basin size (largest basin size divided by the total sum of basin sizes).
    """
    return max(basin_sizes) * 1.0 / sum(basin_sizes)


def get_entropy_of_basin_size_distribution(basin_sizes):
    """
    Compute the Shannon entropy of the basin size distribution.

    This function calculates the Shannon entropy of a probability distribution derived from the basin sizes.
    First, the basin sizes are normalized to form a probability distribution, and then the entropy is computed
    using the formula: H = - sum(p_i * log(p_i)), where p_i is the proportion of the basin size i.

    Parameters:
        basin_sizes (list or array-like): A list where each element represents the size of a basin,
                                           i.e., the number of initial conditions that converge to a particular attractor.

    Returns:
        float: The Shannon entropy of the basin size distribution.
    """
    total = sum(basin_sizes)
    probabilities = [size * 1.0 / total for size in basin_sizes]
    return sum([-np.log(p) * p for p in probabilities])


def get_robustness_from_attractor_dict_exact(attractor_dict, N, n_attractors, left_side_of_truth_table):
    """
    Compute the robustness (coherence) of a Boolean network based on its attractor dictionary.

    This function computes the proportion of neighbors in the Boolean hypercube that, following a synchronous update,
    transition to the same attractor. For each state in the fully sampled state space (left_side_of_truth_table),
    it examines all N neighbors (each obtained by flipping one bit) and counts how many have the same attractor
    as the current state. The robustness is then given as the fraction of such edges over the total number of edges
    in the hypercube (which is 2^(N-1) * N).

    Parameters:
        attractor_dict (dict): A dictionary mapping each state (in decimal representation) to its attractor index.
                                 This dictionary must be computed from a fully sampled state space.
        N (int): The number of nodes in the network.
        n_attractors (int): The total number of attractors (not directly used in computation, but provided for context).
        left_side_of_truth_table (list or array): The full truth table of the network states, where each entry is a numpy array
                                                    representing one state (of length N).

    Returns:
        float: The robustness measure (i.e., the proportion of neighboring states that transition to the same attractor).
    """
    if n_attractors == 1:
        return 1
    b_for_bin2dec = np.array([2**i for i in range(N)])[::-1]
    count_of_neighbors_who_transition_to_same_attractor = 0
    for xdec, x in enumerate(left_side_of_truth_table):
        for i in range(N):
            if x[i] == 0:
                ydec = xdec + b_for_bin2dec[i]
            else:
                continue
            if attractor_dict[xdec] == attractor_dict[ydec]:
                count_of_neighbors_who_transition_to_same_attractor += 1
    return count_of_neighbors_who_transition_to_same_attractor / (2**(N-1) * N)


def get_robustness_measures_and_attractors(F, I, number_different_IC=500):
    """
    Approximate global robustness measures and attractors.

    This function samples the attractor landscape by simulating the network from a number of different initial
    conditions. It computes:
      1. The coherence: the proportion of neighboring states (in the Boolean hypercube) that, after synchronous
         update, transition to the same attractor.
      2. The fragility: a measure of how much the attractor state changes (assumed under synchronous update) in response
         to perturbations.
      3. The final time-step Hamming distance between perturbed trajectories.
    
    In addition, it collects several details about each attractor (such as basin sizes, coherence of each basin, etc.).

    Parameters:
        F (list): List of Boolean functions (truth tables) for each node.
        I (list): List of lists, where I[i] contains the indices of the regulators for node i.
        number_different_IC (int, optional): Number of different initial conditions to sample (default is 500).

    Returns:
        tuple: A tuple containing:
            - attractors (list): List of attractors (each attractor is represented as a list of state decimal numbers).
            - lower_bound_number_of_attractors (int): The lower bound on the number of attractors found.
            - approximate_basin_sizes (list): List of basin sizes for each attractor.
            - approximate_phenotypical_robustness (float): The approximated global robustness (coherence) measure.
            - approximate_fragility (float): The approximated fragility measure (averaged per node).
            - final_hamming_distance_approximation (float): The approximated final Hamming distance measure.
            - approximate_basin_robustness (list): robustness/coherence of each basin.
            - approximate_basin_fragility (list): fragility of each basin.
    """
    import math
    from collections import defaultdict
    def lcm(a, b):
        return abs(a*b) // math.gcd(a, b)
    
    N = len(F)
    
    dictF = dict()
    attractors = []
    ICs_per_attractor_state = []
    basin_sizes = []
    attractor_dict = dict()
    attractor_state_dict = []
    distance_from_attractor_state_dict = []
    counter_phase_shifts = []
    
    height = []
    degrees = list(map(len, I))
    
    b_for_bin2decs = [np.array([2**i for i in range(NN)])[::-1] for NN in range(max(degrees)+1)]
    b_for_bin2dec = np.array([2**i for i in range(N)])[::-1]
    
    robustness_approximation = 0
    fragility_sum = 0
    basin_robustness = defaultdict(float)
    basin_fragility = defaultdict(float)
    final_hamming_distance_approximation = 0
    mean_states_attractors = []
    states_attractors = []
    
    BREAK = False
    for i in range(number_different_IC):
        index_attractors = []
        index_of_state_within_attractor_reached = []
        distance_from_attractor = []
        for j in range(2):
            if j == 0:
                x = np.random.randint(2, size=N)
                xbin = np.dot(x, b_for_bin2dec)
                x_old = x.copy()
            else:
                x = x_old
                random_flipped_bit = np.random.choice(N)
                x[random_flipped_bit] = 1 - x[random_flipped_bit]
                xbin = np.dot(x, b_for_bin2dec)
            queue = [xbin]
            try:
                index_attr = attractor_dict[xbin]
            except KeyError:
                while True:
                    try:
                        fxbin = dictF[xbin]
                    except KeyError:
                        fx = []
                        for jj in range(N):
                            fx.append(F[jj][np.dot(x[I[jj]], b_for_bin2decs[degrees[jj]])])
                        fxbin = np.dot(fx, b_for_bin2dec)
                        dictF.update({xbin: fxbin})
                    try:
                        index_attr = attractor_dict[fxbin]
                        dummy_index_within_attractor_reached = attractor_state_dict[index_attr][fxbin]
                        dummy_distance_from_attractor = distance_from_attractor_state_dict[index_attr][fxbin]
                        attractor_dict.update(list(zip(queue, [index_attr]*len(queue))))
                        attractor_state_dict[index_attr].update(list(zip(queue, [dummy_index_within_attractor_reached]*len(queue))))
                        distance_from_attractor_state_dict[index_attr].update(
                            list(zip(queue, list(range(len(queue) + dummy_distance_from_attractor, dummy_distance_from_attractor, -1))))
                        )
                        break
                    except KeyError:
                        try:
                            index = queue.index(fxbin)
                            index_attr = len(attractors)
                            attractor_dict.update(list(zip(queue, [index_attr]*len(queue))))
                            attractors.append(queue[index:])
                            basin_sizes.append(1)
                            attractor_state_dict.append(dict(zip(queue, [0]*index + list(range(len(attractors[-1])))))
                            )
                            distance_from_attractor_state_dict.append(
                                dict(zip(queue, list(range(index, 0, -1)) + [0]*len(attractors[-1])))
                            )
                            BREAK = False
                            ICs_per_attractor_state.append([0] * len(attractors[-1]))
                            counter_phase_shifts.append([0] * len(attractors[-1]))
                            if BREAK:
                                break
                            if len(attractors[-1]) == 1:
                                states_attractors.append(np.array(dec2bin(queue[index], N)).reshape((1, N)))
                                mean_states_attractors.append(np.array(dec2bin(queue[index], N)))
                            else:
                                states_attractors.append(np.array([dec2bin(state, N) for state in queue[index:]]))
                                mean_states_attractors.append(states_attractors[-1].mean(0))
                            break
                        except ValueError:
                            x = np.array(fx)
                    queue.append(fxbin)
                    xbin = fxbin
                    if BREAK:
                        break
            index_attractors.append(index_attr)
            index_of_state_within_attractor_reached.append(attractor_state_dict[index_attr][xbin])
            distance_from_attractor.append(distance_from_attractor_state_dict[index_attr][xbin])
            basin_sizes[index_attr] += 1
            ICs_per_attractor_state[index_attr][attractor_state_dict[index_attr][xbin]] += 1
        if index_attractors[0] == index_attractors[1]:
            robustness_approximation += 1
            basin_robustness[index_attractors[0]] += 1
            length_phaseshift = max(index_of_state_within_attractor_reached) - min(index_of_state_within_attractor_reached)
            counter_phase_shifts[index_attr][length_phaseshift] += 1
        else:
            fragility_sum += np.sum(np.abs(mean_states_attractors[index_attractors[0]] - mean_states_attractors[index_attractors[1]]))
            basin_fragility[index_attractors[0]] += np.sum(np.abs(mean_states_attractors[index_attractors[0]] - mean_states_attractors[index_attractors[1]]))
            required_n_states = lcm(len(attractors[index_attractors[0]]), len(attractors[index_attractors[1]]))
            index_j0 = index_of_state_within_attractor_reached[0]
            periodic_states_j0 = np.tile(states_attractors[index_attractors[0]], 
                                         (required_n_states // len(attractors[index_attractors[0]]) + 1, 1))[index_j0:(index_j0 + required_n_states), :]
            index_j1 = index_of_state_within_attractor_reached[1]
            periodic_states_j1 = np.tile(states_attractors[index_attractors[1]], 
                                         (required_n_states // len(attractors[index_attractors[1]]) + 1, 1))[index_j1:(index_j1 + required_n_states), :]
            final_hamming_distance_approximation += np.mean(periodic_states_j1 == periodic_states_j0)
            
        height.extend(distance_from_attractor)
    
    lower_bound_number_of_attractors = len(attractors)
    approximate_basin_sizes = basin_sizes
    approximate_phenotypical_robustness = robustness_approximation * 1.0 / number_different_IC
    approximate_fragility = fragility_sum * 1.0 / number_different_IC / N
    
    approximate_basin_robustness = [basin_robustness[index_att] * 2.0 / basin_sizes[index_att]
                                    for index_att in range(len(attractors))]
    approximate_basin_fragility = [basin_fragility[index_att] * 2.0 / basin_sizes[index_att] / N
                                   for index_att in range(len(attractors))]
    
    for index_attr in range(len(attractors)):
        periodic_states_two_periods = np.tile(states_attractors[index_attr], (2, 1))
        for length_phaseshift, num_IC_with_that_phaseshift in enumerate(counter_phase_shifts[index_attr]):
            if num_IC_with_that_phaseshift > 0 and length_phaseshift > 0:
                final_hamming_distance_approximation += num_IC_with_that_phaseshift * np.mean(
                    states_attractors[index_attr] ==
                    periodic_states_two_periods[length_phaseshift:(length_phaseshift + len(attractors[index_attr])), :]
                )
                
    final_hamming_distance_approximation = final_hamming_distance_approximation / number_different_IC
    
    return (attractors, lower_bound_number_of_attractors, approximate_basin_sizes, 
            approximate_phenotypical_robustness, approximate_fragility, final_hamming_distance_approximation,
            approximate_basin_robustness, approximate_basin_fragility)



def adjacency_matrix(I, constants=[], IGNORE_SELFLOOPS=False, IGNORE_CONSTANTS=True):
    """
    Construct the (binary) adjacency matrix from the wiring diagram.

    Given the wiring diagram I (a list of regulator lists for each gene) and a list of constants,
    this function builds an adjacency matrix where each entry m[j, i] is 1 if gene j regulates gene i.
    Self-loops can be optionally ignored, and constant nodes can be excluded.

    Parameters:
        I (list): List of lists, where I[i] contains the indices of regulators for gene i.
        constants (list, optional): List of constant nodes.
        IGNORE_SELFLOOPS (bool, optional): If True, self-loops are ignored.
        IGNORE_CONSTANTS (bool, optional): If True, constant nodes are excluded from the matrix.

    Returns:
        np.array: The binary adjacency matrix.
    """
    n = len(I)
    n_constants = len(constants)
    if IGNORE_CONSTANTS:
        m = np.zeros((n - n_constants, n - n_constants), dtype=int)
        for i in range(len(I)):
            for j in I[i]:
                if j < n - n_constants and (not IGNORE_SELFLOOPS or i != j):
                    m[j, i] = 1
        return m
    else:
        return adjacency_matrix(I, [], IGNORE_CONSTANTS=True)


def get_signed_adjacency_matrix(I, type_of_each_regulation, constants=[], IGNORE_SELFLOOPS=False, IGNORE_CONSTANTS=True):
    """
    Construct the signed adjacency matrix of a Boolean network.

    The signed adjacency matrix assigns +1 for increasing (activating) regulations,
    -1 for decreasing (inhibiting) regulations, and NaN for any other type.

    Parameters:
        I (list): List of lists, where I[i] contains the indices of regulators for gene i.
        type_of_each_regulation (list): List of lists corresponding to the type of regulation ('increasing' or 'decreasing')
                                        for each edge in I.
        constants (list, optional): List of constant nodes.
        IGNORE_SELFLOOPS (bool, optional): If True, self-loops are ignored.
        IGNORE_CONSTANTS (bool, optional): If True, constant nodes are excluded.

    Returns:
        np.array: The signed adjacency matrix.
    """
    n = len(I)
    n_constants = len(constants)
    if IGNORE_CONSTANTS:
        m = np.zeros((n - n_constants, n - n_constants), dtype=int)
        for i, (regulators, type_of_regulation) in enumerate(zip(I, type_of_each_regulation)):
            for j, t in zip(regulators, type_of_regulation):
                if j < n - n_constants and (not IGNORE_SELFLOOPS or i != j):
                    if t == 'increasing':
                        m[j, i] = 1 
                    elif t == 'decreasing':
                        m[j, i] = -1 
                    else:
                        m[j, i] = np.nan
        return m
    else:
        return get_signed_adjacency_matrix(I, type_of_each_regulation, [], IGNORE_CONSTANTS=True)


def get_signed_effective_graph(I, type_of_each_regulation, F, constants=[], IGNORE_SELFLOOPS=False, IGNORE_CONSTANTS=True):
    """
    Construct the signed effective graph of a Boolean network.

    This function computes an effective graph in which each edge is weighted by its effectiveness.
    Effectiveness is obtained via get_edge_effectiveness on the corresponding Boolean function.
    Edges are signed according to the type of regulation ('increasing' or 'decreasing').

    Parameters:
        I (list): List of lists, where I[i] contains the indices of regulators for gene i.
        type_of_each_regulation (list): List of lists specifying the type of regulation for each edge.
        F (list): List of Boolean functions (truth tables) for each node.
        constants (list, optional): List of constant nodes.
        IGNORE_SELFLOOPS (bool, optional): If True, self-loops are ignored.
        IGNORE_CONSTANTS (bool, optional): If True, constant nodes are excluded.

    Returns:
        np.array: The signed effective graph as a matrix of edge effectiveness values.
    """
    n = len(I)
    n_constants = len(constants)
    if IGNORE_CONSTANTS:
        m = np.zeros((n - n_constants, n - n_constants), dtype=float)
        for i, (regulators, type_of_regulation) in enumerate(zip(I, type_of_each_regulation)):
            effectivenesses = get_edge_effectiveness(F[i])
            for j, t, e in zip(regulators, type_of_regulation, effectivenesses):
                if j < n - n_constants and (not IGNORE_SELFLOOPS or i != j):
                    if t == 'increasing':
                        m[j, i] = e
                    elif t == 'decreasing':
                        m[j, i] = -e
                    else:
                        m[j, i] = np.nan
        return m
    else:
        return get_signed_effective_graph(I, type_of_each_regulation, F, [], IGNORE_CONSTANTS=True)


def get_ffls(I, F=None):
    """
    Identify feed-forward loops (FFLs) in a Boolean network and optionally determine their types.

    A feed-forward loop (FFL) is a three-node motif where node i regulates node k both directly and indirectly via node j.
    If F is provided, the function also computes the monotonicity of each regulation in the FFL.

    Parameters:
        I (list): List of lists, where I[i] contains the indices of regulators for gene i.
        F (list, optional): List of Boolean functions (truth tables) for each node.
                             If provided along with F, the types (monotonicities) of the regulations are computed.

    Returns:
        If F is None:
            list: A list of FFLs, each represented as a list [i, k, j] (or similar ordering).
        Otherwise:
            tuple: A tuple (ffls, types), where ffls is a list of FFLs and types is a list of corresponding monotonicity types.
    """
    ffls = []
    types = []
    for i in range(len(I)):
        for j in range(i + 1, len(I)):
            for k in range(len(I)):
                if i == k or j == k:
                    continue
                # Check if there is an FFL: i regulates k and j regulates both i and k.
                if i in I[k] and i in I[j] and j in I[k]:
                    ffls.append([i, j, k])
                    if F is not None:
                        # Compute types if F is provided.
                        # (This example assumes a helper function is_monotonic exists and that I is ordered.)
                        monotonic_i = is_monotonic(F[i], True)[1]
                        monotonic_j = is_monotonic(F[j], True)[1]
                        monotonic_k = is_monotonic(F[k], True)[1]
                        direct = monotonic_k[I[k].index(i)]
                        indirect1 = monotonic_j[I[j].index(i)]
                        indirect2 = monotonic_k[I[k].index(j)]
                        types.append([direct, indirect1, indirect2])
    if F is not None:
        return (ffls, types)
    else:
        return ffls


def get_ffls_from_I(I, types_I=None):
    """
    Identify feed-forward loops (FFLs) in a Boolean network based solely on the wiring diagram.

    The function uses the inverted wiring diagram to identify common targets and returns the FFLs found.
    If types_I (the type of each regulation) is provided, it also returns the corresponding regulation types.

    Parameters:
        I (list): List of lists, where I[i] contains the indices of regulators for gene i.
        types_I (list, optional): List of lists specifying the type (e.g., 'increasing' or 'decreasing') for each regulation.

    Returns:
        If types_I is provided:
            tuple: (ffls, types) where ffls is a list of identified FFLs (each as a list [i, j, k]),
                   and types is a list of corresponding regulation type triplets.
        Otherwise:
            list: A list of identified FFLs.
    """
    all_tfs = list(range(len(I)))
    n_tfs = len(all_tfs)
    all_tfs_dict = dict(zip(all_tfs, list(range(n_tfs))))
    I_inv = [[] for _ in all_tfs]
    for target, el in enumerate(I):
        for regulator in el:
            I_inv[all_tfs_dict[regulator]].append(target)
    ffls = []
    types = []
    for i in range(n_tfs):  # master regulators
        for j in range(n_tfs):
            if i == j or all_tfs[j] not in I_inv[i]:
                continue
            common_targets = list(set(I_inv[i]) & set(I_inv[j]))
            for k in common_targets:
                if all_tfs[j] == k or all_tfs[i] == k:
                    continue
                ffls.append([i, j, k])
                if types_I is not None:
                    direct = types_I[k][I[k].index(all_tfs[i])]
                    indirect1 = types_I[all_tfs[j]][I[all_tfs[j]].index(all_tfs[i])]
                    indirect2 = types_I[k][I[k].index(all_tfs[j])]
                    types.append([direct, indirect1, indirect2])
    if types_I is not None:
        return (ffls, types)
    else:
        return ffls


def get_ffl_type_number(types_vector):
    """
    Compute a numeric type for a feed-forward loop (FFL) based on its regulation types.

    For a given FFL, this function converts the types (strings 'increasing' or 'decreasing') into a numeric code.
    If any type is not in the set {'decreasing', 'increasing'} and 'not essential' is not present, -1 is returned;
    if 'not essential' is present, -2 is returned.

    Parameters:
        types_vector (list): List of regulation type strings for the FFL.

    Returns:
        int: A numeric type representing the FFL, or -1/-2 if types are not as expected.
    """
    if not set(types_vector).issubset(set(['decreasing', 'increasing'])):
        return -1 if 'not essential' not in types_vector else -2
    else:
        dummy = np.array([1 if el == 'increasing' else 0 for el in types_vector])
        nr_type = int(np.dot(dummy, 2**np.arange(len(types_vector))))
    return nr_type


def is_ffl_coherent(types_vector):
    """
    Determine whether a feed-forward loop (FFL) is coherent.

    A coherent FFL is defined (in this context) such that the number of 'increasing' regulations in the FFL is odd.
    If the types are not exclusively 'decreasing' or 'increasing', NaN is returned.

    Parameters:
        types_vector (list): List of regulation type strings for the FFL.

    Returns:
        bool or float: True if the FFL is coherent, False otherwise, or NaN if types are ambiguous.
    """
    if not set(types_vector).issubset(set(['decreasing', 'increasing'])):
        return np.nan
    else:
        dummy = np.array([1 if el == 'increasing' else 0 for el in types_vector])
        COHERENT = (sum(dummy) % 2) == 1
    return COHERENT


def generate_networkx_graph(I, constants, variables):
    """
    Generate a NetworkX directed graph from a wiring diagram.

    Nodes are labeled with variable names (from variables) and constant names (from constants). Edges are added
    from each regulator to its target based on the wiring diagram I.

    Parameters:
        I (list): List of lists, where I[i] contains the indices of regulators for gene i.
        constants (list): List of constant names.
        variables (list): List of variable names.

    Returns:
        networkx.DiGraph: The generated directed graph.
    """
    names = list(variables) + list(constants)
    G = nx.DiGraph()
    G.add_nodes_from(names)
    G.add_edges_from([(names[I[i][j]], names[i]) for i in range(len(variables)) for j in range(len(I[i]))])
    return G


def generate_networkx_graph_from_edges(I, n_variables):
    """
    Generate a NetworkX directed graph from an edge list derived from the wiring diagram.

    Only edges among the first n_variables (excluding constant self-loops) are included.

    Parameters:
        I (list): List of lists, where I[i] contains the indices of regulators for gene i.
        n_variables (int): Number of variable nodes (constants are excluded).

    Returns:
        networkx.DiGraph: The generated directed graph.
    """
    edges = []
    for j, regulators in enumerate(I):
        if j >= n_variables:  # Exclude constant self-loops
            break
        for i in regulators:
            edges.append((i, j))
    return nx.DiGraph(edges)


def simple_cycles(G, max_len=4):
    """
    Generate simple cycles in a directed graph using a variant of Johnson's algorithm.

    This function finds simple cycles (elementary circuits) in graph G with a maximum length of max_len.
    It first yields self-cycles (if any), then iterates through the strongly connected components of G,
    recursively unblocking nodes to yield cycles.

    Parameters:
        G (networkx.DiGraph): A directed graph.
        max_len (int, optional): Maximum length of cycles to consider (default is 4).

    Yields:
        list: A list of nodes representing a simple cycle.
    """
    from collections import defaultdict

    def _unblock(thisnode, blocked, B):
        stack = set([thisnode])
        while stack:
            node = stack.pop()
            if node in blocked:
                blocked.remove(node)
                stack.update(B[node])
                B[node].clear()

    subG = nx.DiGraph(G.edges())
    sccs = [scc for scc in nx.strongly_connected_components(subG) if len(scc) > 1]
    
    # Yield self-cycles and remove them.
    for v in subG:
        if subG.has_edge(v, v):
            yield [v]
            subG.remove_edge(v, v)
    
    while sccs:
        scc = sccs.pop()
        sccG = subG.subgraph(scc)
        startnode = scc.pop()
        path = [startnode]
        len_path = 1
        blocked = set()
        closed = set()
        blocked.add(startnode)
        B = defaultdict(set)
        stack = [(startnode, list(sccG[startnode]))]
        while stack:
            thisnode, nbrs = stack[-1]
            if nbrs and len_path <= max_len:
                nextnode = nbrs.pop()
                if nextnode == startnode:
                    yield path[:]
                    closed.update(path)
                elif nextnode not in blocked:
                    path.append(nextnode)
                    len_path += 1
                    stack.append((nextnode, list(sccG[nextnode])))
                    closed.discard(nextnode)
                    blocked.add(nextnode)
                    continue
            if not nbrs or len_path > max_len:
                if thisnode in closed:
                    _unblock(thisnode, blocked, B)
                else:
                    for nbr in sccG[thisnode]:
                        if thisnode not in B[nbr]:
                            B[nbr].add(thisnode)
                stack.pop()
                path.pop()
                len_path -= 1
        H = subG.subgraph(scc)
        sccs.extend(scc for scc in nx.strongly_connected_components(H) if len(scc) > 1)


def get_type_of_loop(loop, F, I):
    """
    Determine the regulation types along a feedback loop.

    For a given loop (a list of node indices), this function returns a list containing
    the type (e.g., 'increasing' or 'decreasing') of each regulation along the loop.
    The loop is assumed to be ordered such that the first node is repeated at the end.

    Parameters:
        loop (list): List of node indices representing the loop.
        F (list): List of Boolean functions (truth tables) for each node.
        I (list): List of regulator indices for each node.

    Returns:
        list: A list of regulation types corresponding to each edge in the loop.
    """
    n = len(loop)
    dummy = loop[:]
    dummy.append(loop[0])
    res = []
    for i in range(n):
        # Assumes is_monotonic returns a tuple with the monotonicity information.
        res.append(is_monotonic(F[dummy[i+1]], True)[1][list(I[dummy[i+1]]).index(dummy[i])])
    return res


def get_loop_type_number(types_vector):
    """
    Compute a numeric code for a loop based on its regulation types.

    For the given list of regulation types in a loop, this function returns the number
    of 'decreasing' regulations. If the types are not a subset of {'decreasing', 'increasing'},
    it returns -1 (or -2 if 'not essential' is present).

    Parameters:
        types_vector (list): List of regulation type strings.

    Returns:
        int: A numeric code representing the loop type.
    """
    if not set(types_vector).issubset(set(['decreasing', 'increasing'])):
        return -1 if 'not essential' not in types_vector else -2
    else:
        nr_type = int(np.sum([1 if el == 'decreasing' else 0 for el in types_vector]))
    return nr_type


def is_pos_loop(types_vector):
    """
    Determine whether a loop is positive based on its regulation types.

    A positive loop is defined such that the total number of 'decreasing' regulations is even.
    If the types vector contains values other than 'decreasing' or 'increasing', NaN is returned.

    Parameters:
        types_vector (list): List of regulation type strings.

    Returns:
        bool or float: True if the loop is positive, False if negative, or NaN if undefined.
    """
    if not set(types_vector).issubset(set(['decreasing', 'increasing'])):
        return np.nan
    else:
        POSITIVE = (np.sum([1 if el == 'decreasing' else 0 for el in types_vector]) % 2) == 0
    return POSITIVE


def generate_networkx_graph(I, constants, variables):
    """
    Generate a NetworkX directed graph from a wiring diagram.

    Nodes are labeled using the provided variable and constant names.
    Edges are added based on the wiring diagram I.

    Parameters:
        I (list): List of lists, where I[i] contains the indices of regulators for gene i.
        constants (list): List of constant names.
        variables (list): List of variable names.

    Returns:
        networkx.DiGraph: The generated directed graph.
    """
    names = list(variables) + list(constants)
    G = nx.DiGraph()
    G.add_nodes_from(names)
    G.add_edges_from([(names[I[i][j]], names[i]) for i in range(len(variables)) for j in range(len(I[i]))])
    return G


def generate_networkx_graph_from_edges(I, n_variables):
    """
    Generate a NetworkX directed graph from an edge list derived from the wiring diagram.

    Only edges among the first n_variables (i.e., non-constant nodes) are included.

    Parameters:
        I (list): List of lists, where I[i] contains the indices of regulators for gene i.
        n_variables (int): Number of variable nodes (constants are excluded).

    Returns:
        networkx.DiGraph: The resulting directed graph.
    """
    edges = []
    for j, regulators in enumerate(I):
        if j >= n_variables:  # Exclude constant self-loops.
            break
        for i in regulators:
            edges.append((i, j))
    return nx.DiGraph(edges)


















