#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 17:28:53 2025

@author: benco, ckadelka

#Version 6
"""

# -------------------------------------------------------- #
#                   SECTION I - Imports                    #
# -------------------------------------------------------- #
import numpy as np
import itertools
import canalizing_function_toolbox_v18 as can
import networkx as nx

# -------------------------------------------------------- #
#               SECTION II - Conversion Utility            #
# -------------------------------------------------------- #
# II.1: Multistate Integer/Vector
def multi2dec(vector, B):
    column = 1
    decimal = 0
    for i, v in enumerate(vector):
        decimal += v * column
        column *= B[i]
    return decimal

def dec2multi(decimal, B):
    string = ''
    for base_degree in B:
        string += str(decimal % base_degree)
        decimal = int(decimal / base_degree)
    return [int(string[i], B[i]) for i in range(len(string))]

# II.2: Boolean/Multistate Vector
def vec_bin2multi(vec_b, B):
    conversion = [ 0 for _ in range(len(B) + 1) ]
    for i, base in enumerate(B):
        conversion[i + 1] = conversion[i] + base - 1
    vec_ms = [ 0 for _ in range(len(B)) ]
    for i in range(len(B)):
        value = 0
        break_flag = False
        for group_rev_idx in range(conversion[i], conversion[i + 1]):
            idx = conversion[i + 1] - group_rev_idx - 1 + conversion[i]
            if vec_b[idx] > 0:
                if break_flag:
                    return None
                value += 1
            else: 
                break_flag = True
        vec_ms[i] = value
    return vec_ms

def vec_multi2bin(vec_ms, B):
    conversion = [ 0 for _ in range(len(B) + 1) ]
    for i, base in enumerate(B):
        conversion[i + 1] = conversion[i] + base - 1
    vec_b = [ 0 for _ in range(sum(B) - len(B))]
    for i in range(len(vec_ms)):
        value = vec_ms[i]
        idx = conversion[i + 1] - 1
        while idx >= conversion[i] and value > 0:
            vec_b[idx] += 1
            idx -= 1
            value -= 1
    return vec_b

# II.3: Boolean/Multistate Integer
def dec_bin2multi(dec_b, B):
    vector_m = vec_bin2multi(can.dec2bin(dec_b, sum(B) - len(B)), B)
    if vector_m is None:
        return None
    return multi2dec(vector_m, B)

def dec_multi2bin(dec_m, B):
    return can.bin2dec(vec_multi2bin(dec2multi(dec_m, B), B))

# -------------------------------------------------------- #
#               SECTION III - Network Update               #
# -------------------------------------------------------- #
# III.1: Utility
def update_ms_single_node(F, states_regulators, B):
    return F[multi2dec(states_regulators, B)]

# III.2: Synchronous
def update_ms_network_synchronously(F, I, B, X):
    if type(X) == list:
        X = np.array(X)
    if type(B) == list:
        B = np.array(B)
    next_state = np.zeros(len(F), int)
    for i in range(len(F)):
        next_state[i] = update_ms_single_node(F[i], X[I[i]], B[I[i]])
    return next_state

def update_ms_network_synchronously_many_times(F, I, B, X, n_steps):
    for i in range(n_steps):
        current_state = update_ms_network_synchronously(F, I, B, X)
    return current_state

# -------------------------------------------------------- #
#                   SECTION IV - Attractors                #
# -------------------------------------------------------- #
def get_ms_attractors_synchronous(F, I, B, n_sim = 500, initial_sample_points = [], n_steps_timeout = 1000000000000, INITIAL_SAMPLE_POINTS_AS_MUTLISTATE_VECTORS = True):
    func_dict = dict()
    attractors = []
    basin_sizes = []
    attractor_dict = dict()
    state_space = []
    
    sampled_points = []
    n_timeout = 0
    
    node_count = len(F)
    
    INITIAL_SAMPLE_POINTS_EMPTY = ((isinstance(initial_sample_points, np.ndarray) and initial_sample_points.size == 0) or initial_sample_points == [])
    if not INITIAL_SAMPLE_POINTS_EMPTY:
        n_sim = len(initial_sample_points)
    
    for sim_num in range(n_sim):
        if INITIAL_SAMPLE_POINTS_EMPTY:
            sample_state = np.random.randint(B, size = node_count)
            sample_state_decimal = multi2dec(sample_state, B)
            sampled_points.append(sample_state_decimal)
        else:
            if INITIAL_SAMPLE_POINTS_AS_MUTLISTATE_VECTORS:
                sample_state = initial_sample_points[sim_num]
                sample_state_decimal = multi2dec(sample_state)
            else:
                sample_state_decimal = initial_sample_points[sim_num]
                sample_state = np.array(dec2multi(sample_state_decimal, B))
        queue = [sample_state_decimal]
        count = 0
        while count < n_steps_timeout:
            try:
                next_state_decimal = func_dict[sample_state_decimal]
            except KeyError:
                next_state_vector = update_ms_network_synchronously(F, I, B, sample_state)
                next_state_decimal = multi2dec(next_state_vector, B)
                func_dict.update({sample_state_decimal : next_state_decimal})
                sample_state = next_state_vector
            if count == 0:
                state_space.append(next_state_decimal)
            try:
                index_attractor = attractor_dict[next_state_decimal]
                basin_sizes[index_attractor] += 1
                attractor_dict.update(list(zip(queue, [index_attractor] * len(queue))))
                break
            except KeyError:
                try:
                    index = queue.index(next_state_decimal)
                    attractor_dict.update(list(zip(queue[index:], [len(attractors)] * (len(queue) - index))))
                    attractors.append(queue[index:])
                    basin_sizes.append(1)
                    break
                except ValueError:
                    pass
            queue.append(next_state_decimal)
            sample_state_decimal = next_state_decimal
            count += 1
            if count == n_steps_timeout:
                n_timeout += 1
    
    return (attractors, len(attractors), basin_sizes, attractor_dict,
            sampled_points if INITIAL_SAMPLE_POINTS_EMPTY else initial_sample_points,
            state_space, n_timeout)

def get_ms_attractors_synchronous_exact(F, I, B, left_side_of_truth_table = None):
    if type(B) == list:
        B = np.array(B)
    
    N = len(F)
    
    state_size = np.prod(B)
    
    if left_side_of_truth_table is None:
        left_side_of_truth_table = np.array(list(map(np.array, list(dec2multi(d, B) for d in range(state_size)))))
    
    state_space = np.zeros((state_size, N), int)
    for i in range(N):
        count = np.prod(B[I[i]])
        for j in range(count):
            x = np.array(dec2multi(j, B[I[i]]), int)
            if F[i][j]:
                state_space[np.all(left_side_of_truth_table[:, I[i]] == x, axis = 1), i] = F[i][j]
    
    dictF = dict()
    for i in range(state_size):
        dictF[i] = multi2dec(update_ms_network_synchronously(F, I, B, dec2multi(i, B)), B)
    
    attractors = []
    basin_sizes = []
    attractor_dict = dict()
    for state_decimal in range(state_size):
        queue = [ state_decimal ]
        while True:
            next_state_decimal = dictF[state_decimal]
            try:
                index_attr = attractor_dict[next_state_decimal]
                basin_sizes[index_attr] += 1
                attractor_dict.update(list(zip(queue, [index_attr] * len(queue))))
                break
            except KeyError:
                try:
                    index = queue.index(next_state_decimal)
                    attractor_dict.update(list(zip(queue, [len(attractors)] * len(queue))))
                    attractors.append(queue[index:])
                    basin_sizes.append(1)
                    break
                except ValueError:
                    pass
            queue.append(next_state_decimal)
            state_decimal = next_state_decimal
    return (attractors, len(attractors), basin_sizes, attractor_dict, state_space)

def get_bms_attractors_synchronous(F, I, B, nsim=500, initial_sample_points=[], n_steps_timeout=1000000000000,
                               INITIAL_SAMPLE_POINTS_AS_BINARY_VECTORS=True):
    dictF = dict()
    attractors = []
    basin_sizes = []
    attr_dict = dict()
    state_space = []
    
    sampled_points = []
    n_timeout = 0
    
    N = len(F)
    
    INITIAL_SAMPLE_POINTS_EMPTY = ((isinstance(initial_sample_points, np.ndarray) and initial_sample_points.size == 0) or initial_sample_points == [])
    if not INITIAL_SAMPLE_POINTS_EMPTY:
        nsim = len(initial_sample_points)
    
    for i in range(nsim):
        if INITIAL_SAMPLE_POINTS_EMPTY:
            x = np.random.randint(2, size=N)
            while vec_bin2multi(x, B) == None: # prevent sampling invalid states
                x = np.random.randint(2, size=N)
            xdec = can.bin2dec(x)
            sampled_points.append(xdec)
        else:
            if INITIAL_SAMPLE_POINTS_AS_BINARY_VECTORS:
                x = initial_sample_points[i]
                xdec = can.bin2dec(x)
            else:
                xdec = initial_sample_points[i]
                x = np.array(can.dec2bin(xdec, N))
        queue = [xdec]
        count = 0
        while count < n_steps_timeout:
            try:
                fxdec = dictF[xdec]
            except KeyError:
                fx = can.update_network_synchronously(F, I, x)
                fxdec = can.bin2dec(fx)
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

def get_bms_attractors_synchronous_exact(F, I, B, left_side_of_truth_table=None):
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
        if dec_bin2multi(xdec, B) == None:
            continue
        
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

# -------------------------------------------------------- #
#                   SECTION V - Analysis                   #
# -------------------------------------------------------- #
# V.1: Booleanization
def booleanize(F, I, B):
    # Create a helper conversion array, which stores the starting index for every
    # group of Boolean nodes that cumulatively represent the same multistate node
    conversion = [ 0 for _ in range(len(B) + 1) ]
    for i, base in enumerate(B):
        conversion[i + 1] = conversion[i] + base - 1
    
    # Default value to -1, so it will be very obvious if something goes wrong as we are
    # expecting I_bool to be a matrix (array of arrays) of non-negative integers, so
    # if I_bool directly contains a -1 then something went wrong here
    I_bool = [ -1 for _ in range(conversion[len(conversion) - 1]) ]
    
    for i, reg_idxs in enumerate(I):
        # We need to generate the interior array, _I
        _I = []
        for multi_idx in reg_idxs:
            for bool_idx in range(conversion[multi_idx], conversion[multi_idx + 1]):
                _I.append(bool_idx)
        
        # We ultimately want _I to be a numpy array, not a list
        _I = np.array(_I, int)
        
        # Every boolean node within the consecutive pair of indices in the conversion array
        # is representative of the same multistate node, such that for some integer A:
        #       i.e. a pair (A, A + 1) is a binary node in the multistate network
        #       i.e. a pair (A, A + 2) is a ternary node in the multistate network
        # As they are representative of the same mulitstate node, we can simply assign
        # them the same regulatory nodes, as the logic will never extend beyond these
        # confines, even if a specific node may not require all of these regulators itself
        for j in range(conversion[i], conversion[i + 1]):
            I_bool[j] = _I
    
    # Default value to -1, so it will be very obvious if something goes wrong as we are
    # expecting F_bool to be a matrix (array of arrays) of 0s and 1s, so if F_bool directly
    # contains a -1 then something went wrong here
    F_bool = [ -1 for _ in range(len(I_bool)) ]
    
    for i_b in range(len(I_bool)):
        # The Boolean indices will increase at a faster rate than the multistate indices
        # as we have defined a multistate node with N states to be represented by N-1 Boolean
        # nodes. As such, we need to determine our position in the multistate network (i_m) from
        # our position in the Booleanized network (i_b).
        i_m = 0
        while i_b >= conversion[i_m + 1]:
            i_m += 1
        
        # We need to generate the interior array, _F, which has a length of 2^(# of regulators)
        _F = [ 0 for _ in range(2 ** len(I_bool[i_b])) ]
        
        sub_base = B[I[i_m]]
        for i_f_m in range(len(F[i_m])):
            _F[dec_multi2bin(i_f_m, sub_base)] = int(F[i_m][i_f_m] >= conversion[i_m + 1] - i_b)
        F_bool[i_b] = np.array(_F, int)
    
    return (F_bool, I_bool)

def unbooleanize(F, I, B):
    # Guarantee the debooleanization is valid
    if sum(B) - len(B) != len(I):
        print(f"Error: Unexpected groupings for multistate bases (expected multistate bases for {len(I)} nodes, received for {sum(B) - len(B)} nodes)")
        return
    
    idxm = 0
    idxb = 0
    conversion = dict()
    for b in B:
        for x in range(idxb, idxb + b - 1):
            conversion[x] = idxm
        idxb += b - 1
        idxm += 1
    
    I_multi = [ set() for _ in range(len(B)) ]
    for i, reg_idxs in enumerate(I):
        _I = set()
        for bool_idx in reg_idxs:
            _I.add(conversion[bool_idx])
        I_multi[conversion[i]] = I_multi[conversion[i]].union(_I)
    for i in range(len(I_multi)):
        if len(I_multi[i]) == 0:
            print("Error: Unexpected regulator matrix, received empty row")
            return
        I_multi[i] = np.array(list(I_multi[i]), int)
    
    F_multi = [ 0 for _ in range(len(I_multi)) ]
    for i_b in range(len(I)):
        i_m = conversion[i_b]
        
        _F = [ 0 for _ in range(np.prod(B[I_multi[i_m]], 0, np.int64)) ]
        
        i_f_m = 0
        B_Im_c_ib = B[I_multi[conversion[i_b]]]
        for i_f_b in range(len(F[i_b])):
            d_b2m = dec_bin2multi(i_f_b, B_Im_c_ib)
            if d_b2m == None:
                continue
            _F[d_b2m] = F[i_b][i_f_b]
            i_f_m += 1
        
        F_multi[i_m] += np.array(_F, int)
    
    return (F_multi, I_multi)

# V.2: State Space
def compare_state_space_MSN_BMSN(F, I, B):
    ss_multi = set(multi2dec(state, B) for state in get_ms_attractors_synchronous_exact(F, I, B)[4])
    
    Fb, Ib = booleanize(F, I, B)
    ss_bin = set(multi2dec(vec_bin2multi(state, B), B) for state in can.get_attractors_synchronous_exact(Fb, Ib)[4])
    
    expected = []
    unexpected = []
    found = []
    
    for state in ss_multi:
        if state in ss_bin:
            found.append(state)
            ss_bin.remove(state)
        else:
            expected.append(state)
    for state in ss_bin:
        unexpected.append(state)
    
    return (found, expected, unexpected)

# V.3: Attractors
def compare_attractors_basin_sizes(attr_A, basin_A, attr_B, basin_B, key_separator_char = '→'):
    # Define helper method to convert list of cyclically ordered integers into
    # unique string representation regardless of phase shift in ordering
    def __attr2key__(attr):
        init = sorted(attr)[0]
        bfor = ""
        aftr = ""
        flag = False
        for nmbr in attr:
            if nmbr == init:
                flag = True
                continue
            if flag:
                aftr += str(nmbr) + key_separator_char
            else:
                bfor += str(nmbr) + key_separator_char
        return (str(init) + key_separator_char + aftr + bfor)
    
    attr_dict = dict()
    
    attr_idx_A = dict()
    attr_idx_B = dict()
    
    basin_sizes_A = dict()
    basin_sizes_B = dict()
    
    # Get dictionary of all attractors within superset and
    # respective indices in subsets A and B
    for i, attr in enumerate(attr_A):
        key = __attr2key__(attr)
        attr_dict[key] = True
        attr_idx_A[key] = i
    for i, attr in enumerate(attr_B):
        key = __attr2key__(attr)
        attr_dict[key] = True
        attr_idx_B[key] = i
    
    for key in attr_dict.keys():
        try:
            basin_sizes_A[key] = basin_A[attr_idx_A[key]]
        except KeyError:
            basin_sizes_A[key] = 0
        try:
            basin_sizes_B[key] = basin_B[attr_idx_B[key]]
        except KeyError:
            basin_sizes_B[key] = 0
    
    return (basin_sizes_A, basin_sizes_B)

def compare_attractors_MSN_BMSN(F, I, B, num_sample_points = 500, key_separator_char = '→'):
    attr_m, _, basin_m, _, _, _, _ = get_ms_attractors_synchronous(F, I, B, n_sim = num_sample_points)
    
    Fb, Ib = booleanize(F, I, B)
    attr_b, _, basin_b, _, _, _, _ = can.get_attractors_synchronous(Fb, Ib, nsim = num_sample_points)
    for i, attr in enumerate(attr_b):
        attr_b[i] = [ dec_bin2multi(state, B) for state in attr ]
    
    return compare_attractors_basin_sizes(attr_m, basin_m, attr_b, basin_b, key_separator_char)

def compare_attractors_MSN_BMSN_exact(F, I, B, key_separator_char = '→'):
    attr_m, _, basin_m, _, _ = get_ms_attractors_synchronous_exact(F, I, B)
    
    Fb, Ib = booleanize(F, I, B)
    attr_b, _, basin_b, _, _ = can.get_attractors_synchronous_exact(Fb, Ib)
    for i, attr in enumerate(attr_b):
        attr_b[i] = [ dec_bin2multi(state, B) for state in attr ]
    
    return compare_attractors_basin_sizes(attr_m, basin_m, attr_b, basin_b, key_separator_char)

# -------------------------------------------------------- #
#           SECTION VI - Random Network Generation         #
# -------------------------------------------------------- #
# VI.1: Utility
def f_from_ms_expression(expr, reg_bases):
    expr_split = expr.replace('(', ' ( ').replace(')', ' ) ').split(' ')
    var = []
    var_dict = dict()
    num_var = 0
    for i, el in enumerate(expr_split):
        if el not in [ '', ' ', '(', ')', 'and', 'or', 'not', 'AND', 'OR', 'NOT', '&', '|', '~', '+', '-', '*', '%', '>', '>=', '==', '<=', '<' ] and not el.isdigit():
            try:
                new_var = var_dict[el]
            except KeyError:
                new_var = 'x[%i]' % num_var
                var_dict[el] = new_var
                var.append(el)
                num_var += 1
            expr_split[i] = new_var
        elif el in ['AND', 'OR', 'NOT']:
            expr_split[i] - el.lower()
    expr = ' '.join(expr_split)
    f = []
    f_size = 1
    for b in reg_bases:
        f_size *= b
    for x in range(f_size):
        x = dec2multi(x, reg_bases)
        f.append(int(eval(expr)))
    return f, var

# VI.2: Generation
def random_non_degenerated_ms_function(in_degree, my_base, reg_bases):
    while True:
        f = np.array(np.floor(np.random.random(np.prod(reg_bases)) * my_base), int)
        #if not is_ms_degenerated(f, reg_bases):
        return f

def random_linear_ms_function(in_degree, my_base, reg_bases):
    expr = '(%s) %% %i' % (' + '.join(['(%i + x%i)' % (np.random.randint(0, my_base), y) for y in range(in_degree)]), my_base)
    return f_from_ms_expression(expr, reg_bases)[0]

def random_MSN(N, in_degree = 2, base = 3, STRONGLY_CONNECTED = True, in_degree_distribution = 'const',
               UNIFORM_BASE_DISTRIBUTION = False, list_x = [], NO_SELF_REGULATION = True, edges_wiring_diagram = None, LINEAR = False):
    if in_degree_distribution in [ 'constant', 'const', 'dirac', 'delta' ]:
        if type(in_degree) in [ list, np.array ]:
            try:
                assert np.all([ type(el) in [ int, np.int64 ] for el in in_degree ])
                assert len(in_degree) == N
                assert min(in_degree) >= 1
                assert max(in_degree) <= N
                ns = np.array(in_degree[:], int)
            except AssertionError:
                print("Error: A vector was submitted for the in-degree.\nTo use a user-defined in-degree vector, ensure that it is an N-dimensional vector where each element is an integer between 1 and the N.")
                return
        else:
            try:
                assert type(in_degree) in [ int, np.int64 ]
                assert in_degree >= 1
                assert in_degree <= N
                ns = np.ones(N, int) * in_degree
            except AssertionError:
                print("Error: The in-degree must be a single integer (or N-dimensional vector of integers) between 1 and N when using a constant degree distribution.")
                return
    elif in_degree_distribution == 'uniform':
        if type(in_degree) in [ list, np.array ]:
            try:
                assert np.all([ type(el) in [ int, np.int64 ] for el in in_degree ])
                assert len(in_degree) == N
                assert min(in_degree) >= 1
                assert max(in_degree) <= N
                ns = np.array(in_degree[:])
            except AssertionError:
                print("Error: A vector was submitted for the in-degree.\nEnsure that you are providing an N-dimensional vector where each element is an integer between 1 and N representing the upper bound of a uniform degree distribution (lower bound == 1).")
                return
        else:
            try:
                assert type(in_degree) in [ int, np.int64 ]
                assert in_degree >= 1
                assert in_degree <= N
                ns = np.ones(N, int) * in_degree
            except AssertionError:
                print("Error: The in-degree must be a single integer (or N-dimensional vector of integers) between 1 and N representing the upper bound of a uniform degree distribution (lower bound == 1).")
                return
    elif in_degree_distribution == 'poisson':
        if type(in_degree) in [ list, np.array ]:
            try:
                assert np.all([ type(el) in [ int, np.int64, float, np.float64 ] for el in in_degree ])
                assert len(in_degree) == N
                assert min(in_degree) >= 1
                assert max(in_degree) <= N
                ns = np.array(in_degree[:])
            except AssertionError:
                print("Error: A vector was submitted for the in-degree.\nEnsure that the in-degree is an N-dimensional vector where each element is > 0 and represents the Poisson parameter.")
                return
        else:
            try:
                assert type(in_degree) in [ int, np.int64, float, np.float64 ]
                assert in_degree >= 1
                assert in_degree <= N
                ns = np.ones(N, int) * in_degree
            except AssertionError:
                print("Error: The in-degree must be a single number (or N-dimensional vector) > 0 representing the Poisson parameter.")
                return
    else:
        print("Error: None of the predefined in-degree distrbutions were chosen.\nTo use a user-defined in-degree vector, use the input 'in_degree' to submit an N-dimensional vector where each element must be between 1 and N.")
        return
    
    if edges_wiring_diagram is None:
        while True: # continue generation until a strongly connected graph is generated
            if in_degree_distribution == 'uniform':
                ns = 1 + np.random.randint(in_degree - 1, None, N, np.int64)
            elif in_degree_distribution == 'poisson':
                ns = np.random.poisson(in_degree, N)
                ns[ns == 0] = 1
                ns[ns > N - int(NO_SELF_REGULATION)] = N - int(NO_SELF_REGULATION)
            edges_wiring_diagram = can.random_edge_list(N, ns, NO_SELF_REGULATION)
            if STRONGLY_CONNECTED:
                G = nx.from_edgelist(edges_wiring_diagram, create_using = nx.MultiDiGraph())
                if not nx.is_strongly_connected(G):
                    continue
            break
    else:
        try:
            assert len(set(np.array(edges_wiring_diagram).flatten())) == N
        except AssertionError:
            print("Number of nodes provided in edges_wiring_diagram != N")
            return
        ns = np.zeros(N, int)
        for target in np.array(edges_wiring_diagram)[:, 1]:
            ns[target] += 1
    
    if UNIFORM_BASE_DISTRIBUTION:
        if type(base) in [ list, np.array ]:
            try:
                assert np.all([ type(el) in [ int, np.int64 ] for el in base ])
                assert len(base) == N
                assert min(base) >= 2
                B = np.array([ 2 + np.random.randint(base[i] - 1, None, None) for i in range(N) ], int)
            except AssertionError:
                print("Error: A vector was submitted for the base.\nEnsure that you are providing an N-dimensional vector where each element is an integer >= 2 representing the upper bound of a uniform distribution (lower bound == 2).")
                return
        else:
            try:
                assert type(base) in [ int, np.int64 ]
                assert base >= 2
                B = 2 + np.random.randint(base - 1, None, N, np.int64)
            except AssertionError:
                print("Error: The base must be a single integer (or N-dimensional vector of integers) >= 2 representing the upper bound of a uniform distribution (lower bound == 2).")
                return
    else:
        if type(base) in [ list, np.array ]:
            try:
                assert np.all([ type(el) in [ int, np.int64 ] for el in base ])
                assert len(base) == N
                assert min(base) >= 2
                B = np.array(base[:], int)
            except AssertionError:
                print("Error: A vector was submitted for the base.\nTo use a user-defined base vector, ensure that it is an N-dimensional vector where each element is an integer >= 2.")
                return
        else:
            try:
                assert type(base) in [ int, np.int64 ]
                assert base >= 2
                B = np.ones(N, int) * base
            except AssertionError:
                print("Error: The base must be a single integer (or N-dimensional vector of integers) >= 2 when using a constant distribution.")
                return
    
    I = [ [] for _ in range(N) ]
    for edge in edges_wiring_diagram:
        I[edge[1]].append(edge[0])
    for i in range(N):
        I[i] = np.sort(I[i])
    
    F = []
    for i in range(N):
        if LINEAR:
            F.append(random_linear_ms_function(ns[i], B[i], B[I[i]]))
        else:
            F.append(random_non_degenerated_ms_function(ns[i], B[i], B[I[i]]))
    
    return F, I, B, ns

# -------------------------------------------------------- #
#                SECTION VII - Benchmarking                #
# -------------------------------------------------------- #
import time

def _calc_attr(F, I, B, num_sample_points = 1000, key_separator_char = '→'):
    s_m = time.process_time_ns()
    attr_m, _, basin_m, _, _, _, _ = get_ms_attractors_synchronous(F, I, B, n_sim = num_sample_points)
    e_m = time.process_time_ns()
    
    s_b1 = time.process_time_ns()
    Fb, Ib = booleanize(F, I, B)
    e_b1 = time.process_time_ns()
    
    s_b2 = time.process_time_ns()
    attr_b, _, basin_b, _, _, _, _ = can.get_attractors_synchronous(Fb, Ib, nsim = num_sample_points)
    e_b2 = time.process_time_ns()
    
    s_b3 = time.process_time_ns()
    for i, attr in enumerate(attr_b):
        attr_b[i] = [ dec_bin2multi(state, B) for state in attr ]
    e_b3 = time.process_time_ns()
    
    return((e_m - s_m), (e_b1 - s_b1), (e_b2 - s_b2), (e_b3 - s_b3))

def _calc_attr_ex(F, I, B, key_separator_char = '→'):
    s_m = time.process_time_ns()
    attr_m, _, basin_m, _, _ = get_ms_attractors_synchronous_exact(F, I, B)
    e_m = time.process_time_ns()
    
    s_b1 = time.process_time_ns()
    Fb, Ib = booleanize(F, I, B)
    e_b1 = time.process_time_ns()
    
    s_b2 = time.process_time_ns()
    attr_b, _, basin_b, _, _ = can.get_attractors_synchronous_exact(Fb, Ib)
    e_b2 = time.process_time_ns()
    
    s_b3 = time.process_time_ns()
    for i, attr in enumerate(attr_b):
        attr_b[i] = [ dec_bin2multi(state, B) for state in attr ]
    e_b3 = time.process_time_ns()
    
    return((e_m - s_m), (e_b1 - s_b1), (e_b2 - s_b2), (e_b3 - s_b3))

def _trial(N, in_degree, base, force_strongly_connected, n_sim = -1):
    F, I, B, _ = random_MSN(N, in_degree, base, STRONGLY_CONNECTED=force_strongly_connected, UNIFORM_BASE_DISTRIBUTION=False, in_degree_distribution='const')
    
    if n_sim <= 0:
        time_m, time_booleanize, time_b, time_debooleanize = _calc_attr_ex(F, I, B)
    else:
        time_m, time_booleanize, time_b, time_debooleanize = _calc_attr(F, I, B, n_sim)
    results = { "Multistate" : time_m, "BooleanizeFunction" : time_booleanize, "Booleanized" : time_b, "DebooleanizeResult" : time_debooleanize, "BooleanizedTotal" : (time_booleanize + time_b + time_debooleanize) }
    return results