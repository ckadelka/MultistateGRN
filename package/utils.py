#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# """
# Utility functions used throughout BoolForge.

# The :mod:`~boolforge.utils` module includes low-level operations for binary and
# decimal conversions, truth table manipulations, and combinatorial helper
# functions. These utilities are used internally by
# :class:`~boolforge.BooleanFunction` and :class:`~boolforge.BooleanNetwork`
# classes to enable efficient representation and analysis of Boolean functions
# and networks.

# Notes
# -----
# Most functions in this module are intended for internal use and are not part of
# the stable public API.

# Examples
# --------
# >>> import boolforge
# >>> boolforge.bin2dec([1, 0, 1])
# 5
# >>> boolforge.dec2bin(5, 3)
# array([1, 0, 1])
# """

import random as _py_random
from collections.abc import Sequence
import numpy as np
from numpy.random import Generator as _NPGen, RandomState as _NPRandomState, SeedSequence, default_rng
#from typing import Tuple

_LOGIC_MAP = {
    "AND": "&",
    "and": "&",
    "&&": "&",
    "&": "&",
    "OR": "|",
    "or": "|",
    "||": "|",
    "|": "|",
    "NOT": "~",
    "not": "~",
    "!": "~",
    "~": "~",
}

_COMPARE_OPS = {"==", "!=", ">=", "<=", ">", "<"}

_ARITH_OPS = {"+", "-", "*", "%"}

# def _require_cana():
#     try:
#         import cana.boolean_node
#         return cana.boolean_node
#     except ModuleNotFoundError as e:
#         raise ImportError(
#             "This functionality requires CANA. "
#             "Install it with `pip install cana`."
#         ) from e
        

def _coerce_rng(
    rng : int | _NPGen | _NPRandomState | _py_random.Random | None = None
) -> _NPGen:
    """
    Coerce a variety of RNG-like inputs to a NumPy ``Generator``.

    Parameters
    ----------
    rng : int | np.random.Generator | np.random.RandomState | random.Random | None, optional
        Random number generator or seed specification.

        - ``None``: return ``np.random.default_rng()``.
        - ``int``: interpreted as a seed for ``default_rng``.
        - ``np.random.Generator``: returned unchanged.
        - ``np.random.RandomState``: converted via ``SeedSequence``.
        - ``random.Random``: converted via ``SeedSequence``.

    Returns
    -------
    np.random.Generator
        A NumPy random number generator.

    Raises
    ------
    TypeError
        If ``rng`` is not one of the supported types.

    Notes
    -----
    This function provides a unified RNG interface across BoolForge by
    normalizing legacy and standard-library RNGs to the modern NumPy
    ``Generator`` API.

    Conversion from ``RandomState`` and ``random.Random`` is performed by
    extracting entropy and initializing a ``SeedSequence``. This preserves
    reproducibility while avoiding direct reliance on deprecated RNG APIs.
    """
    if rng is None:
        return default_rng()
    if isinstance(rng, _NPGen):
        return rng
    if isinstance(rng, (int, np.integer)):
        return default_rng(int(rng))
    if isinstance(rng, _NPRandomState):
        # derive robust entropy from the legacy RNG
        entropy = rng.randint(0, 2**32, size=4, dtype=np.uint32)
        return default_rng(SeedSequence(entropy))
    if isinstance(rng, _py_random.Random):
        entropy = [_py_random.Random(rng.random()).getrandbits(32) for _ in range(4)]
        # simpler: entropy = [rng.getrandbits(32) for _ in range(4)]
        return default_rng(SeedSequence(entropy))
    raise TypeError(f"Unsupported rng type: {type(rng)!r}")


def _is_number(token: str) -> bool:
    """Return True if token is a pure numeric literal."""
    try:
        float(token)
        return True
    except ValueError:
        return False
    

# def is_float(element: object) -> bool:
#     """
#     Check whether an object can be coerced to a float.

#     Parameters
#     ----------
#     element : object
#         Object to test for float coercibility.

#     Returns
#     -------
#     bool
#         True if ``element`` can be converted to ``float`` without raising an
#         exception, False otherwise.

#     Notes
#     -----
#     This function tests coercibility, not type membership. For example,
#     numeric strings and integers return True.

#     Examples
#     --------
#     >>> is_float(3)
#     True
#     >>> is_float(3.14)
#     True
#     >>> is_float("2.7")
#     True
#     >>> is_float("abc")
#     False
#     >>> is_float(None)
#     False
#     """
#     try:
#         float(element)
#         return True
#     except (TypeError, ValueError):
#         return False

def mix2dec(vector : Sequence[int], radices : Sequence[int]) -> int:
    """
    Convert a mixed-radix vector to an integer.
    
    Parameters
    ----------
    vector : Sequence of int
        Digits ordered from least significant digit to most significant digit.
    radices : Sequence of int
        Radices for each digit in the vector.
    
    Returns
    -------
    int
        Integer represented by the mixed-radix vector.
    """
    column = 1
    decimal = 0
    for i, v in enumerate(vector):
        decimal += v * column
        column *= radices[i]
    return decimal

def dec2mix(decimal : int, radices : Sequence[int]) -> list[int]:
    """
    Convert a nonnegative integer into a mixed-radix vector.

    Parameters
    ----------
    decimal : int
        Nonnegative integer to convert.
    radices : Sequence[int]
        Radices for each digit in the resulting vector.

    Returns
    -------
    list[int]
        Digits ordered from least significant digit to most significant digit.
    """
    vector = []
    for radix in radices:
        vector.append(decimal % radix)
        decimal //= radix
    return vector

def mix2bin(vector : Sequence[int], radices : Sequence[int]) -> list[int]:
    """
    Convert a mixed-radix vector to a binary vector.
    
    Parameters
    ----------
    vector : Sequence of int
        Digits ordered from least significant digit to most significant digit.
    radices : Sequence of int
        Radices for each digit in the vector.
    
    Returns
    -------
    list of int
        Binary digits (0 or 1), ordered from most significant bit to least
        significant bit.
    """
    decimal = mix2dec(vector, radices)
    binstr = bin(decimal)[2:].zfill(sum(radices) - len(radices))
    return [ int(bit) for bit in binstr ]

def bin2mix(binary_vector : Sequence[int], radices : Sequence[int]) -> list[int]:
    """
    Convert a binary vector to a mixed-radix vector.
    
    Parameters
    ----------
    binary_vector : Sequence of int
        Binary digits (0 or 1), ordered from most significant bit to least
        significant bit.
    radices : Sequence of int
        Radices for each digit in the resulting mixed-radix vector.
    
    Returns
    -------
    list of int
        Digits ordered from least significant digit to most significant digit.
    """
    decimal = 0
    for bit in binary_vector:
        decimal = (decimal << 1) | bit
    return dec2mix(decimal, radices)


left_side_of_truth_tables = {}

def get_left_side_of_truth_table(R : Sequence[int]) -> np.ndarray:
    """
    Return the left-hand side of a Boolean truth table.

    The left-hand side is the binary representation of all input
    combinations for ``N`` Boolean variables, ordered lexicographically.

    Parameters
    ----------
    R : Sequence[int]

    Returns
    -------
    np.ndarray
        Array of shape ``(2**N, N)`` with entries in ``{0, 1}``. Columns are
        ordered from most significant bit to least significant bit.
    """
    R = np.prod(R)
    if R in left_side_of_truth_tables:
        left_side_of_truth_table = left_side_of_truth_tables[R]
    else:
        left_side_of_truth_table = np.arange(R, dtype=np.uint64)[:, None]
        left_side_of_truth_tables[R] = left_side_of_truth_table
    return left_side_of_truth_table


def find_all_indices(arr: list, el: object) -> list[int]:
    """
    Find all indices of a given element in a sequence.

    Parameters
    ----------
    arr : list
        Sequence to search.
    el : object
        Element to locate.

    Returns
    -------
    list of int
        Indices ``i`` such that ``arr[i] == el``.

    Raises
    ------
    ValueError
        If ``el`` does not occur in ``arr``.

    Examples
    --------
    >>> find_all_indices([1, 2, 1, 3], 1)
    [0, 2]
    >>> find_all_indices(['a', 'b', 'a'], 'a')
    [0, 2]
    """
    res: list[int] = []
    for i, a in enumerate(arr):
        if a == el:
            res.append(i)

    if not res:
        raise ValueError("Element not found in sequence")

    return res


# def check_if_empty(my_list: list | np.ndarray) -> bool:
#     """
#     Check whether a list or NumPy array is empty.

#     Parameters
#     ----------
#     my_list : list or np.ndarray
#         Sequence to check.

#     Returns
#     -------
#     bool
#         True if ``my_list`` is empty, False otherwise.

#     Notes
#     -----
#     For NumPy arrays, emptiness is determined by ``size == 0``.
#     For Python lists, emptiness is determined by equality to ``[]``.
#     """
#     if isinstance(my_list, np.ndarray):
#         return my_list.size == 0
#     return my_list == []
    
    
# def is_list_or_array_of_ints(
#     x: list | np.ndarray,
#     required_length: int | None = None,
# ) -> bool:
#     """
#     Check whether a list or NumPy array contains only integers.

#     Parameters
#     ----------
#     x : list or np.ndarray
#         Sequence to check.
#     required_length : int or None, optional
#         If provided, require that ``x`` has exactly this length.

#     Returns
#     -------
#     bool
#         True if ``x`` is a list of ``int`` / ``np.integer`` or a NumPy array
#         with integer dtype, and (if specified) has length ``required_length``.
#         False otherwise.

#     Notes
#     -----
#     - For Python lists, each element is checked individually.
#     - For NumPy arrays, the dtype is checked using ``np.issubdtype``.
#     - One-dimensional arrays are required when ``required_length`` is given.
#     """
#     # Case 1: Python list
#     if isinstance(x, list):
#         return (
#             (required_length is None or len(x) == required_length)
#             and all(isinstance(el, (int, np.integer)) for el in x)
#         )

#     # Case 2: NumPy array
#     if isinstance(x, np.ndarray):
#         return (
#             (required_length is None or x.shape == (required_length,))
#             and np.issubdtype(x.dtype, np.integer)
#         )

#     return False


# def is_list_or_array_of_floats(
#     x: list | np.ndarray,
#     required_length: int | None = None,
# ) -> bool:
#     """
#     Check whether a list or NumPy array contains only floating-point numbers.

#     Parameters
#     ----------
#     x : list or np.ndarray
#         Sequence to check.
#     required_length : int or None, optional
#         If provided, require that ``x`` has exactly this length.

#     Returns
#     -------
#     bool
#         True if ``x`` is a list of ``float`` / ``np.floating`` or a NumPy array
#         with floating-point dtype, and (if specified) has length
#         ``required_length``. False otherwise.

#     Notes
#     -----
#     - For Python lists, each element is checked individually.
#     - For NumPy arrays, the dtype is checked using ``np.issubdtype``.
#     - One-dimensional arrays are required when ``required_length`` is given.
#     """
#     # Case 1: Python list
#     if isinstance(x, list):
#         return (
#             (required_length is None or len(x) == required_length)
#             and all(isinstance(el, (float, np.floating)) for el in x)
#         )

#     # Case 2: NumPy array
#     if isinstance(x, np.ndarray):
#         return (
#             (required_length is None or x.shape == (required_length,))
#             and np.issubdtype(x.dtype, np.floating)
#         )

#     return False


# def bool_to_poly(
#     f: list,
#     variables: list[str] | None = None,
#     prefix: str = '',
# ) -> str:
#     """
#     Convert a Boolean function from truth-table form to disjunctive normal form.

#     The returned expression is a non-reduced disjunctive normal form (DNF),
#     expressed as a sum of monomials corresponding to truth-table entries where
#     the function evaluates to 1.

#     Parameters
#     ----------
#     f : list
#         Boolean function values ordered according to the standard truth-table
#         convention. The length of ``f`` must be ``2**n`` for some integer ``n``.
#     variables : list of str or None, optional
#         Variable names to use in the expression. If None or if the length does
#         not match the required number of variables, default names
#         ``['x0', 'x1', ..., 'x{n-1}']`` are used.
#     prefix : str, optional
#         Prefix for automatically generated variable names. Ignored if
#         ``variables`` is provided with the correct length.

#     Returns
#     -------
#     str
#         Boolean expression in non-reduced disjunctive normal form. Returns
#         ``'0'`` if the function is identically zero.

#     Notes
#     -----
#     - Variables are ordered from most significant bit to least significant bit.
#     - Each monomial corresponds to a single truth-table row where ``f == 1``.
#     - No simplification or reduction of the DNF is performed.

#     Examples
#     --------
#     >>> bool_to_poly([0, 1, 1, 0])
#     '(1 - x0) * x1 + x0 * (1 - x1)'
#     """
#     len_f = len(f)
#     n = int(np.log2(len_f))
#     if variables is None or len(variables) != n:
#         prefix = 'x'
#         variables = [prefix + str(i) for i in range(n)]

#     left_side_of_truth_table = get_left_side_of_truth_table(n)
#     num_values = 2 ** n
#     text = []

#     for i in range(num_values):
#         if f[i] == True:
#             monomial = ' * '.join(
#                 [
#                     v if entry == 1 else f'(1 - {v})'
#                     for v, entry in zip(variables, left_side_of_truth_table[i])
#                 ]
#             )
#             text.append(monomial)

#     if text:
#         return ' + '.join(text)
#     return '0'


def f_from_expression(expr, R, max_degree = 16):
    expr = expr.replace('(', ' ( ').replace(')', ' ) ')
    raw_tokens = expr.split()
    tokens = []
    variables = []
    seen = set()
    
    for token in raw_tokens:
        if token in {"(", ")"}:
            tokens.append(token)
            continue
        if token in _LOGIC_MAP:
            tokens.append(_LOGIC_MAP[token])
            continue
        if token in _COMPARE_OPS:
            tokens.append(token)
            continue
        if token in _ARITH_OPS:
            tokens.append(token)
            continue
        if _is_number(token):
            tokens.append(token)
            continue
        if token not in seen:
            seen.add(token)
            variables.append(token)
        token.appnd(token)
    n = len(variables)
    if n > max_degree:
        return np.array([], np.uint8), np.array(variables)
    safe_map = { var: f"v{i}" for i, var in enumerate(variables) }
    safe_tokens = [
        safe_map[token] if token in safe_map else token for token in tokens
    ]
    expr_mod = " ".join(safe_tokens)
    truth_table = get_left_side_of_truth_table(n, R)
    local_dict = {
        safe_map[var]: truth_table[:, i].astype(np.int64) for i, var in enumerate(variables)
    }
    try:
        result = eval(expr_mod, {"__builtins__":None},local_dict)
    except Exception as e:
        raise ValueError(f"Error evaluating expression: \n{expr}\nParsed as:\n{expr_mod}\nError: {e}")
    result = np.asarray(result)
    if n == 0:
        result = np.array([int(result)], np.int64)
    else:
        result = result.astype(np.int64)
    result %= np.prod(R)
    return result.astype(np.uint8), np.array(variables)


# def flatten(l: Sequence[Sequence[object]]) -> list[object]:
#     """
#     Flatten a sequence of sequences by one level.

#     Parameters
#     ----------
#     l : list or np.ndarray
#         Sequence whose elements are themselves iterable.

#     Returns
#     -------
#     list
#         A flat list containing the elements of each sub-sequence in ``l``,
#         in order.

#     Notes
#     -----
#     This function performs a single-level flattening only. Nested sequences
#     deeper than one level are not recursively flattened.

#     Examples
#     --------
#     >>> flatten([[1, 2], [3, 4]])
#     [1, 2, 3, 4]
#     >>> flatten(np.array([[1, 2], [3, 4]]))
#     [1, 2, 3, 4]
#     """
#     return [item for sublist in l for item in sublist]


# def hamming_weight_to_ncf_layer_structure(
#     n: int,
#     w: int,
# ) -> list[int]:
#     """
#     Compute the canalizing layer structure of a nested canalizing function (NCF)
#     from its Hamming weight.

#     For nested canalizing functions, there is a bijection between the (odd)
#     Hamming weight ``w`` and the canalizing layer structure, with ``w`` and
#     ``2**n - w`` corresponding to the same structure.

#     Parameters
#     ----------
#     n : int
#         Number of input variables of the NCF.
#     w : int
#         Odd Hamming weight of the NCF.

#     Returns
#     -------
#     list of int
#         Canalizing layer structure ``[k_1, ..., k_r]``.

#     Raises
#     ------
#     TypeError
#         If ``w`` is not an integer.
#     ValueError
#         If ``w`` is outside ``[1, 2**n - 1]`` or if ``w`` is even.

#     Notes
#     -----
#     - All nested canalizing functions have odd Hamming weight.
#     - The binary expansion of ``w`` (with ``n`` bits) determines the layer
#       structure.
      
#     References
#     ----------
#     Kadelka, C., Kuipers, J., & Laubenbacher, R. (2017). 
#     The influence of canalization on the robustness of Boolean networks. 
#     Physica D: Nonlinear Phenomena, 353, 39-47.
#     """
#     if not isinstance(w, (int, np.integer)):
#         raise TypeError("Hamming weight w must be an integer")

#     if not (1 <= w <= 2**n - 1):
#         raise ValueError("Hamming weight w must satisfy 1 <= w <= 2**n - 1")

#     if w % 2 == 0:
#         raise ValueError("Hamming weight w must be odd for nested canalizing functions")

#     if w == 1:
#         return [n]

#     w_bin = dec2bin(w, n)

#     current_el = w_bin[0]
#     layer_structure_NCF = [1]

#     for el in w_bin[1:-1]:
#         if el == current_el:
#             layer_structure_NCF[-1] += 1
#         else:
#             layer_structure_NCF.append(1)
#             current_el = el

#     layer_structure_NCF[-1] += 1
#     return layer_structure_NCF

# # ===================== #
# #   Modular BoolForge   #
# # ===================== #

# import math
# import networkx as nx

# def merge_state_representation(x : int | Sequence[int], 
#                                y : int | Sequence[int],
#                                num_nodes : int | Sequence[int]
#     ) -> int | Sequence[int]:
#     """
#     Combine two state representations into a single decimal representation.
    
#     Parameters
#     ----------
#     x : int or sequence of int
#         First state. Can be a single integer or a pair of integers.
    
#     y : int or sequence of int
#         Second state. Can be a single integer or a pair of integers.
    
#     b : int or sequence of int
#         Bit size of y. Must match the structure of y (int or pair of ints).
    
#     Returns
#     -------
#     result : int or tuple of int
#         Combined state representation. Returns an int if both x and y
#         are integers. Returns a tuple of two ints if either x or y is a tuple/list.
#     """

#     is_pair_x = isinstance(x, Sequence)
#     is_pair_y = isinstance(y, Sequence)
#     if is_pair_x:
#         if is_pair_y:
#             return ((x[0] << num_nodes[0]) | y[0], (x[1] << num_nodes[1]) | y[1]) 
#         return (x[0], (x[1] << num_nodes) | y)
#     elif is_pair_y:
#         return (y[0], (x << num_nodes[1]) | y[1])
#     return (x << num_nodes) | y

# def get_product_of_attractors(attrs_1 : Sequence[Sequence[int | Sequence[int]]],
#     attrs_2 : Sequence[Sequence[int | Sequence[int]]],
#     bits : int | Sequence[int]) -> list:
#     """
#     Compute the product of two sets of attractors by combining their states.
    
#     Parameters
#     ----------
#     attrs_1 : sequence of sequences of int or sequence of sequences of pairs of ints
#         First set of attractors. Each attractor is a list of states.
    
#     attrs_2 : sequence of sequences of int or sequence of sequences of pairs of ints
#         Second set of attractors. Each attractor is a list of states.
    
#     bits : int or pair of ints
#         Bit size of states in attrs_2. Used when merging states.
    
#     Returns
#     -------
#     result : sequence of sequences of int or sequence of sequences of pairs of ints
#         Product set of attractors obtained by merging each attractor
#         from attrs_1 with each attractor from attrs_2.
#     """

#     attractors = []
#     for attr1 in attrs_1:
#         attr = []
#         for attr2 in attrs_2:
#             m = len(attr1)
#             n = len(attr2)
#             for i in range(math.lcm(*[m, n])):
#                 attr.append(merge_state_representation(attr1[i % m], attr2[i % n], bits))
#         attractors.append(attr)
#     return attractors

# def compress_trajectories(trajectories : tuple[Sequence[int], int],
#     num_nodes : int) -> nx.DiGraph:
#     """
#     Compress multiple trajectories into a single directed graph.
    
#     Each trajectory is represented by a prefix (non-periodic states)
#     and a cycle (periodic states). Nodes are merged when identical
#     prefixes or cycles occur across trajectories.
    
#     Parameters
#     ----------
#     trajectories : tuple of (sequence of int, int)
#         List of trajectories. Each trajectory is a tuple containing
#         a list of decimal states and the length of its periodic cycle.
    
#     num_nodes : int
#         Number of nodes in the network. Used to format node labels as binary strings.
    
#     Returns
#     -------
#     G : networkx.DiGraph
#         Directed graph representing all merged trajectories.
#     """

#     # Helper method: determine the 'canon' ordering of a periodic pattern.
#     # The canon ordering is the phase such that the lowest states come first
#     # without changing the relative ordering of the states.
#     def _canon_cycle_(pattern):
#         return min([ tuple(pattern[i:] + pattern[:i]) for i in range(len(pattern)) ])
    
#     # Helper method: determine which offset a given pattern is from the canon
#     # ordering. That is, how much the pattern has been phased relative to the
#     # canon ordering.
#     def _cycle_offset_(pattern, canon):
#         pattern = list(pattern)
#         canon = list(canon)
#         len_pattern = len(pattern)
#         for offset in range(len_pattern):
#             if canon[offset:] + canon[:offset] == pattern:
#                 return offset
#         raise ValueError("Pattern does not match canonical rotations")
    
#     G = nx.DiGraph()
#     next_id = 0
#     cycle_nodes = {}
#     prefix_merge = {}
#     for states, period in trajectories:
#         len_traj = len(states)
#         # First look through the non-periodic component of the trajectory,
#         # also referred to in this code as the 'prefix' of the trajectory
#         len_pref = len_traj - period
#         pref_ids = []
#         for i in range(len_pref):
#             # Determine if this prefix can be merged elsewhere into the graph
#             future = states[i:]
#             prefix_tail = future[:-period]
#             pattern = future[-period:]
#             canon = _canon_cycle_(pattern)
#             entry_offset = _cycle_offset_(pattern, canon)
#             signature = (tuple(prefix_tail), canon, entry_offset)
#             # If so, merge the it and mark the node as initial
#             if signature in prefix_merge:
#                 node_id = prefix_merge[signature]
#                 if i == 0:
#                     G.nodes[node_id]["StIn"] = True
#             # Otherwise, make a new initial node
#             else:
#                 node_id = next_id
#                 prefix_merge[signature] = node_id
#                 G.add_node(next_id, StIn=(i == 0),
#                     NLbl=(str(dec2bin(states[i], num_nodes)).replace(' ', '').replace(',', '').replace('[', '').replace(']', '')))
#                 pref_ids.append(next_id)
#                 next_id += 1
#             pref_ids.append(node_id)
#         # Once prefix nodes are added, create edges
#         for i in range(len(pref_ids) - 1):
#             if pref_ids[i] != pref_ids[i+1]:
#                 G.add_edge(pref_ids[i], pref_ids[i+1])
#         # Second look through the periodic component of the trajectory,
#         # also referred to in this code as the 'cycle' of the trajectory
#         cycle = states[-period:]
#         key = _canon_cycle_(cycle)
#         # If we have found a new cycle, add it to the graph
#         if key not in cycle_nodes:
#             ids = []
#             for s in key:
#                 # Create nodes based off of the canon ordering to ensure
#                 # predictable ordering in case we need to reference
#                 # this cycle again for another trajectory
#                 G.add_node(next_id, StIn=False,
#                     NLbl=(str(dec2bin(s, num_nodes)).replace(' ', '').replace(',', '').replace('[', '').replace(']', '')))
#                 ids.append(next_id)
#                 next_id += 1
#             # Once nodes are added, add in edges
#             for a, b in zip(ids, ids[1:]):
#                 G.add_edge(a, b)
#             G.add_edge(ids[-1], ids[0])
#             cycle_nodes[key] = ids
#         # For a trajectory without a prefix, mark the first state of the trajectory
#         # within the cycle as an initial node
#         if len_pref == 0:
#             G.nodes()[cycle_nodes[key][_cycle_offset_(cycle, key)]]["StIn"] = True
#         # Otherwise, we need to add an edge between the prefix and cycle
#         else:
#             G.add_edge(pref_ids[-1], cycle_nodes[key][_cycle_offset_(cycle, key)])
#     return G

# def product_of_trajectories(compressed_trajectory_graph_1 : nx.DiGraph,
#     compressed_trajectory_graph_2 : nx.DiGraph) -> nx.DiGraph:
#     """
#     Compute the product of two compressed trajectory graphs, following the
#     premise of equal reachability.
    
#     The resulting graph contains all combinations of nodes from
#     the two input graphs, with edges representing all possible
#     successor pairs.
    
#     Parameters
#     ----------
#     compressed_trajectory_graph_1 : networkx.DiGraph
#         First compressed trajectory graph.
    
#     compressed_trajectory_graph_2 : networkx.DiGraph
#         Second compressed trajectory graph.
    
#     Returns
#     -------
#     G : networkx.DiGraph
#         Directed graph representing the product of the two input graphs.
#     """

#     _initial_1 = []
#     _initial_2 = []
#     for n in compressed_trajectory_graph_1.nodes:
#         if compressed_trajectory_graph_1.nodes[n]["StIn"]:
#             _initial_1.append(n)
#     for n in compressed_trajectory_graph_2.nodes:
#         if compressed_trajectory_graph_2.nodes[n]["StIn"]:
#             _initial_2.append(n)
#     G = nx.DiGraph()
#     starting = []
#     for n1 in _initial_1:
#         for n2 in _initial_2:
#             starting.append((n1, n2))
#             G.add_node((n1, n2), StIn=compressed_trajectory_graph_1.nodes[n1]["StIn"] and compressed_trajectory_graph_2.nodes[n2]["StIn"],
#                 NLbl=f"{compressed_trajectory_graph_1.nodes[n1]['NLbl']}{compressed_trajectory_graph_2.nodes[n2]['NLbl']}")
#     stack = starting[:]
#     visited = set(starting)
#     while stack:
#         u1, u2 = stack.pop()
#         for v1 in compressed_trajectory_graph_1.successors(u1):
#             for v2 in compressed_trajectory_graph_2.successors(u2):
#                 new_pair = (v1, v2)
#                 if new_pair not in G:
#                     G.add_node(new_pair, StIn=False,
#                         NLbl=f"{compressed_trajectory_graph_1.nodes[v1]['NLbl']}{compressed_trajectory_graph_2.nodes[v2]['NLbl']}")
#                 G.add_edge((u1, u2), new_pair)
#                 if new_pair not in visited:
#                     visited.add(new_pair)
#                     stack.append(new_pair)
#     return G

# def plot_trajectory(compressed_trajectory_graph : nx.DiGraph,
#     show : bool = True):
#     """
#     Visualize a compressed trajectory graph using a layered layout.
    
#     Initial states are highlighted with a box. Layers are computed
#     based on weakly connected components to improve readability.
    
#     Parameters
#     ----------
#     compressed_trajectory_graph : networkx.DiGraph
#         Directed graph of compressed trajectories.
    
#     show : bool, default=True
#         Whether to call ``plt.show()`` at the end.
#     """
#     import matplotlib.pyplot as plt
    
#     def layout_tree(G, root, x0, y0, dx, pos, visited):
#         children = [c for c in G.predecessors(root) if c not in visited]
#         if not children:
#             return
#         width = dx * (len(children) - 1)
#         xs = [x0 - width/2 + i*dx for i in range(len(children))]
        
#         for child, x in zip(children, xs):
#             pos[child] = (x, y0 - 1)
#             visited.add(child)
#             layout_tree(G, child, x, y0 - 1, dx/1.5, pos, visited)
#         return
    
#     G = compressed_trajectory_graph.copy()
    
#     components = list(nx.weakly_connected_components(G))
#     fig, axes = plt.subplots(
#         nrows=len(components),
#         figsize=(10, 5 * len(components)),
#         squeeze=False
#     )

#     labels = nx.get_node_attributes(G, "NLbl")
#     initial = nx.get_node_attributes(G, "StIn")

#     for idx, comp in enumerate(components):
#         ax = axes[idx][0]
#         sub_nodes = list(comp)
#         SG = G.subgraph(sub_nodes).copy()
#         pos = {}

#         # Find a cycle in the component
#         start = sub_nodes[0]
#         visited = {}
#         v = start
#         while v not in visited:
#             visited[v] = True
#             succ = list(SG.successors(v))
#             if not succ:
#                 break
#             v = succ[0]

#         # Build the cycle
#         cycle = [v]
#         succ = list(SG.successors(v))
#         if succ:
#             u = succ[0]
#             while u != v:
#                 cycle.append(u)
#                 u = next(iter(SG.successors(u)))

#         # Place cycle nodes in a circle
#         n = len(cycle)
#         for i, node in enumerate(cycle):
#             angle = 2 * math.pi * i / n
#             pos[node] = (2 * np.cos(angle), 2 * np.sin(angle))

#         # Layout trees hanging off the cycle
#         visited = set(cycle)
#         for node in cycle:
#             layout_tree(SG, node, pos[node][0], pos[node][1], dx=1.5, pos=pos, visited=visited)
        
#         nx.draw_networkx(
#             SG,
#             pos,
#             ax=ax,
#             node_size=1200,
#             node_color="#00000000",
#             arrows=True,
#             arrowstyle="-|>",
#             arrowsize=8
#         )
        
#         ax.invert_yaxis() # flip, so the graph points downward
        
#         normal = {}
#         boxed = {}
#         for n in SG.nodes():
#             if initial.get(n, False):
#                 boxed[n] = labels[n]
#             else:
#                 normal[n] = labels[n]
        
#         nx.draw_networkx_labels(
#             SG, pos, labels=normal, font_size=12,
#             bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="white", lw=1),
#             ax=ax
#         )
#         nx.draw_networkx_labels(
#             SG, pos, labels=boxed, font_size=12,
#             bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="black", lw=1),
#             ax=ax
#         )
        
#         ax.axis("equal")
#         ax.axis("off")

#     plt.tight_layout()
#     if show:
#         plt.show()
    
#     return fig