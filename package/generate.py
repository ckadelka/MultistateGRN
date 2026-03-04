#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module provides functions for generating random Boolean functions and
Boolean networks with specified structural and dynamical properties.

The :mod:`~boolforge.generate` module enables the systematic creation of
Boolean functions and networks that satisfy particular constraints, such as
specified canalization depth, sensitivity range, bias, or connectivity.
Generated instances can be used for statistical analysis, benchmarking, or
simulation studies.

Several generation routines leverage Numba acceleration for efficient sampling
and evaluation of large function spaces. While Numba is **recommended** to
achieve near-native performance, it is **not required** for functionality; all
functions have pure Python fallbacks.

This module complements :mod:`~boolforge.boolean_function` and
:mod:`~boolforge.boolean_network` by facilitating reproducible generation of
synthetic test cases and large ensembles of random networks.

Example
-------
>>> import boolforge
>>> boolforge.random_function(n=3)
>>> boolforge.random_network(N=5, n=2)
"""


##Imports
import math

import numpy as np
import networkx as nx
from collections.abc import Sequence

from .boolean_function import BooleanFunction
from .boolean_network import BooleanNetwork, WiringDiagram
from . import utils

__all__ = [
    "random_function",
    "random_function_with_bias",
    "random_function_with_exact_hamming_weight",
    "random_degenerate_function",
    "random_non_degenerate_function",
    "random_non_canalizing_function",
    "random_non_canalizing_non_degenerate_function",
    "random_parity_function",
    "random_k_canalizing_function",
    "random_k_canalizing_function_with_specific_layer_structure",
    "random_NCF",
    "random_network",
    "random_null_model",
    "random_wiring_diagram",
    "random_edge_list",
    "random_degrees",
    "rewire_wiring_diagram",
]


## Random function generation
def _validate_bias(bias: float) -> None:
    if not isinstance(bias, (float, int, np.floating)):
        raise TypeError("bias must be a float")
    if not (0.0 <= bias <= 1.0):
        raise ValueError("bias must be in [0, 1]")

def _validate_absolute_bias(absolute_bias: float) -> None:
    if not isinstance(absolute_bias, (float, int, np.floating)):
        raise TypeError("absolute_bias must be a float")
    if not (0.0 <= absolute_bias <= 1.0):
        raise ValueError("absolute_bias must be in [0, 1]")

def _validate_hamming_weight(
    n: int,
    hamming_weight: int,
    *,
    exact_depth: bool,
) -> None:
    if not isinstance(hamming_weight, (int, np.integer)):
        raise TypeError("hamming_weight must be an integer")
    if not (0 <= hamming_weight <= 2**n):
        raise ValueError("hamming_weight must satisfy 0 <= hamming_weight <= 2**n")

    if exact_depth and not (1 < hamming_weight < 2**n - 1):
        raise ValueError(
            "If exact_depth=True and depth=0, hamming_weight must be in "
            "{2, 3, ..., 2**n - 2}. "
            "Functions with weights 0, 1, 2**n-1, 2**n are canalizing."
        )



def random_function(
    n: int,
    depth: int = 0,
    exact_depth: bool = False,
    uniform_over_functions: bool = True,
    layer_structure: list[int] | None = None,
    parity: bool = False,
    allow_degenerate_functions: bool = False,
    bias: float = 0.5,
    absolute_bias: float = 0,
    use_absolute_bias: bool = False,
    hamming_weight: int | None = None,
    *,
    rng = None,
) -> BooleanFunction:
    """
    Generate a random Boolean function under flexible structural constraints.

    This function acts as a high-level generator that unifies several common
    ensembles of Boolean functions, including parity functions, canalizing
    functions of specified depth or layer structure, functions with fixed
    Hamming weight, and biased random functions. The first applicable
    generation rule (in the order described below) is applied.

    Selection logic (first applicable rule is used)

    1. If ``parity`` is True, return a random parity function
       (see ``random_parity_function``).

    2. Else, if ``layer_structure`` is provided, return a Boolean function
       with the specified canalizing layer structure using
       ``random_k_canalizing_function_with_specific_layer_structure``.
       Exactness of the canalizing depth is controlled by ``exact_depth``.

    3. Else, if ``depth > 0``, return a k-canalizing function with
       ``k = min(depth, n)`` using ``random_k_canalizing_function``.
       If ``exact_depth`` is True, the function has exactly this depth;
       otherwise, its canalizing depth is at least ``k``.

       If ``uniform_over_functions`` is True, canalizing layer structures are
       sampled uniformly at random (up to the imposed constraints).
       If False, canalized outputs are sampled independently and uniformly
       as bitstrings, which biases the distribution toward more symmetric
       layer structures.

    4. Else, if ``hamming_weight`` is provided, repeatedly sample Boolean
       functions with the specified Hamming weight until additional
       constraints implied by ``exact_depth`` and
       ``allow_degenerate_functions`` are satisfied.

    5. Else, generate a random Boolean function using a Bernoulli model with
       either:

       - fixed bias ``bias``, or
       - an automatically chosen bias determined by ``absolute_bias`` if
         ``use_absolute_bias`` is True.

       Additional constraints on canalization and degeneracy are enforced
       depending on ``exact_depth`` and ``allow_degenerate_functions``.

    Parameters
    ----------
    n : int
        Number of input variables. Must be a positive integer.
    depth : int, optional
        Requested canalizing depth. Used only if ``layer_structure`` is None
        and ``depth > 0``. If ``exact_depth`` is True, the function has exactly
        this canalizing depth (clipped at ``n``); otherwise, its depth is at
        least ``depth``. Default is 0.
    exact_depth : bool, optional
        Enforce exact canalizing depth where applicable. If ``depth == 0``,
        setting ``exact_depth=True`` enforces that the function is
        non-canalizing. Default is False.
    uniform_over_functions : bool, optional
        If True (default), canalizing layer structures are sampled uniformly
        at random in canalizing-function branches. If False, canalized outputs
        are sampled independently as bitstrings, inducing a bias toward more
        symmetric structures.
    layer_structure : list[int] or None, optional
        Explicit canalizing layer structure ``[k1, ..., kr]``. If provided,
        this takes precedence over ``depth``. Default is None.
    parity : bool, optional
        If True, ignore all other options and return a random parity function.
        Default is False.
    allow_degenerate_functions : bool, optional
        If True, functions with non-essential variables may be returned in
        random-generation branches. If False, non-degenerate functions are
        enforced whenever possible. Default is False.
    bias : float, optional
        Probability of a 1 when sampling truth-table entries independently.
        Used only if ``use_absolute_bias`` is False and no other branch applies.
        Must lie in ``[0, 1]``. Default is 0.5.
    absolute_bias : float, optional
        Absolute deviation from 0.5 used to determine the bias when
        ``use_absolute_bias`` is True. The bias is chosen uniformly from
        ``{0.5*(1 - absolute_bias), 0.5*(1 + absolute_bias)}``. Must lie in
        ``[0, 1]``. Default is 0.
    use_absolute_bias : bool, optional
        If True, ignore ``bias`` and determine the bias using
        ``absolute_bias``. Default is False.
    hamming_weight : int or None, optional
        If provided, enforce that the Boolean function has exactly this many
        ones in its truth table. Additional constraints are enforced depending
        on ``exact_depth`` and ``allow_degenerate_functions``. Default is None.
    rng : int, numpy.random.Generator, numpy.random.RandomState, random.Random, or None, optional
        Random number generator or seed specification. Passed to
        ``utils._coerce_rng``.

    Returns
    -------
    BooleanFunction
        A randomly generated Boolean function of arity ``n``.

    Raises
    ------
    TypeError
        If parameters have invalid types.
    ValueError
        If parameter values or combinations are invalid.

    Notes
    -----
    For any fixed combination of parameters, this function samples **uniformly
    at random** from the set of Boolean functions satisfying the corresponding
    constraints. Non-uniformity arises only when explicitly requested via
    ``uniform_over_functions=False``.

    Extremely biased functions are often degenerate or highly canalizing;
    under restrictive parameter choices, some branches may reject repeatedly
    before returning a valid function.

    Examples
    --------
    >>> # Unbiased, non-degenerate random function
    >>> f = random_function(n=3)

    >>> # Function with canalizing depth at least 2
    >>> f = random_function(n=5, depth=2)

    >>> # Function with exact canalizing depth 2
    >>> f = random_function(n=5, depth=2, exact_depth=True)

    >>> # Function with a specific canalizing layer structure
    >>> f = random_function(n=6, layer_structure=[2, 1])

    >>> # Parity function
    >>> f = random_function(n=4, parity=True)

    >>> # Fixed Hamming weight with non-canalizing and non-degenerate constraints
    >>> f = random_function(
    ...     n=5,
    ...     hamming_weight=10,
    ...     exact_depth=True,
    ...     allow_degenerate_functions=False
    ... )
    """

    if not isinstance(n, (int, np.integer)) or n <= 0:
        raise ValueError("n must be a positive integer")

    rng = utils._coerce_rng(rng)

    # ------------------------------------------------------------
    # Parity branch (highest priority)
    # ------------------------------------------------------------
    if parity:
        return random_parity_function(n, rng=rng)

    # ------------------------------------------------------------
    # Layer structure branch
    # ------------------------------------------------------------
    if layer_structure is not None:
        return random_k_canalizing_function_with_specific_layer_structure(
            n,
            layer_structure,
            exact_depth=exact_depth,
            rng=rng,
        )

    # ------------------------------------------------------------
    # Canalizing depth branch
    # ------------------------------------------------------------
    if not isinstance(depth, (int, np.integer)) or depth < 0:
        raise ValueError("depth must be a nonnegative integer")
    
    if depth > 0:
        return random_k_canalizing_function(
            n,
            min(depth, n),
            exact_depth=exact_depth,
            uniform_over_functions=uniform_over_functions,
            rng=rng,
        )

    # ------------------------------------------------------------
    # Fixed Hamming weight branch
    # ------------------------------------------------------------
    if hamming_weight is not None:
        _validate_hamming_weight(n, hamming_weight, exact_depth=exact_depth)

        while True:
            f = random_function_with_exact_hamming_weight(
                n, hamming_weight, rng=rng
            )

            if allow_degenerate_functions and exact_depth:
                if not f.is_canalizing():
                    return f

            elif allow_degenerate_functions:
                return f

            elif exact_depth:
                if not f.is_canalizing() and not f.is_degenerate():
                    return f

            else:
                if not f.is_degenerate():
                    return f

    # ------------------------------------------------------------
    # Bias-based random generation
    # ------------------------------------------------------------
    if use_absolute_bias:
        _validate_absolute_bias(absolute_bias)
        bias_of_function = rng.choice(
            [0.5 * (1 - absolute_bias), 0.5 * (1 + absolute_bias)]
        )
    else:
        _validate_bias(bias)
        bias_of_function = bias

    if allow_degenerate_functions:
        if exact_depth:
            return random_non_canalizing_function(
                n, bias_of_function, rng=rng
            )
        else:
            return random_function_with_bias(
                n, bias_of_function, rng=rng
            )

    else:
        if exact_depth:
            return random_non_canalizing_non_degenerate_function(
                n, bias_of_function, rng=rng
            )
        else:
            return random_non_degenerate_function(
                n, bias_of_function, rng=rng
            )


def random_function_with_bias(
    n: int,
    bias: float = 0.5,
    *,
    rng=None,
) -> BooleanFunction:
    """
    Generate a random Boolean function with a specified bias.

    The Boolean function is represented by its truth table of length
    ``2**n``, where each entry is independently set to 1 with probability
    ``bias`` and to 0 otherwise.

    Parameters
    ----------
    n : int
        Number of Boolean variables.
    bias : float, optional
        Probability that a given truth-table entry equals 1. Default is 0.5.
    rng : int, np.random.Generator, np.random.RandomState, random.Random, or None, optional
        Random number generator or seed specification. Passed to
        ``utils._coerce_rng``.

    Returns
    -------
    BooleanFunction
        Random Boolean function with the specified bias.
    """
    rng = utils._coerce_rng(rng)
    return BooleanFunction._from_f_unchecked(
        np.array(rng.random(2**n) < bias, dtype=int)
    )


def random_function_with_exact_hamming_weight(
    n: int,
    hamming_weight: int,
    *,
    rng=None,
) -> BooleanFunction:
    """
    Generate a random Boolean function with a fixed Hamming weight.

    The Boolean function is represented by its truth table of length
    ``2**n``, containing exactly ``hamming_weight`` entries equal to 1.
    All such functions are sampled uniformly at random.

    Parameters
    ----------
    n : int
        Number of Boolean variables.
    hamming_weight : int
        Number of truth-table entries equal to 1.
    rng : int, np.random.Generator, np.random.RandomState, random.Random, or None, optional
        Random number generator or seed specification. Passed to
        ``utils._coerce_rng``.

    Returns
    -------
    BooleanFunction
        Random Boolean function with exactly ``hamming_weight`` ones in its
        truth table.

    Raises
    ------
    TypeError
        If ``hamming_weight`` is not an integer.
    ValueError
        If ``hamming_weight`` is not in the range ``[0, 2**n]``.
    """
    rng = utils._coerce_rng(rng)

    if not isinstance(hamming_weight, (int, np.integer)):
        raise TypeError("hamming_weight must be an integer")

    if not (0 <= hamming_weight <= 2**n):
        raise ValueError("hamming_weight must satisfy 0 <= hamming_weight <= 2**n")

    one_indices = rng.choice(2**n, hamming_weight, replace=False)
    f = np.zeros(2**n, dtype=int)
    f[one_indices] = 1

    return BooleanFunction._from_f_unchecked(f)


def random_parity_function(
    n: int,
    *,
    rng=None,
) -> BooleanFunction:
    """
    Generate a random parity Boolean function.

    A parity Boolean function evaluates to the parity (sum modulo 2) of all
    input variables, optionally shifted by a constant. This function returns
    either the parity function or its complement, chosen uniformly at random.

    Parameters
    ----------
    n : int
        Number of Boolean variables.
    rng : int, np.random.Generator, np.random.RandomState, random.Random, or None, optional
        Random number generator or seed specification. Passed to
        ``utils._coerce_rng``.

    Returns
    -------
    BooleanFunction
        Random parity Boolean function on ``n`` variables.

    Raises
    ------
    ValueError
        If ``n`` is not a positive integer.

    Notes
    -----
    - The returned function is either
      ``x1 XOR x2 XOR ... XOR xn`` or its complement.
    - All variables are included symmetrically.
    - Parity functions are never canalizing. All variables must always be known
      to determine the output; they have maximal average sensitivity.

    Examples
    --------
    >>> f = random_parity_function(3)
    >>> sum(f)
    4
    """
    if not isinstance(n, (int, np.integer)) or n <= 0:
        raise ValueError("n must be a positive integer")

    rng = utils._coerce_rng(rng)

    # choose parity or its complement
    val = rng.integers(2)

    f = np.zeros(2**n, dtype=np.uint8)
    for i in range(2**n):
        if i.bit_count() % 2 == val:
            f[i] = 1

    return BooleanFunction._from_f_unchecked(f)



def random_non_degenerate_function(
    n: int,
    bias: float = 0.5,
    *,
    rng=None,
) -> BooleanFunction:
    """
    Generate a random non-degenerate Boolean function.

    A Boolean function is non-degenerate if every variable is essential, i.e.,
    the function depends on all ``n`` input variables. Functions are sampled
    repeatedly from the Bernoulli(bias) ensemble until a non-degenerate
    function is obtained.

    Parameters
    ----------
    n : int
        Number of Boolean variables.
    bias : float, optional
        Probability that a truth-table entry equals 1. Default is 0.5.
    rng : int, np.random.Generator, np.random.RandomState, random.Random, or None, optional
        Random number generator or seed specification. Passed to
        ``utils._coerce_rng``.

    Returns
    -------
    BooleanFunction
        Random non-degenerate Boolean function.

    Raises
    ------
    ValueError
        If ``n`` is not a positive integer.
    ValueError
        If ``bias`` is not strictly between 0 and 1.

    Notes
    -----
    - For moderate bias values, almost all Boolean functions are non-degenerate.
    - Extremely biased functions are very likely to be degenerate, which may
      lead to long rejection-sampling times.
    """
    if not isinstance(n, (int, np.integer)) or n <= 0:
        raise ValueError("n must be a positive integer")

    if not isinstance(bias, (float, np.floating)) or not (0.0 < bias < 1.0):
        raise ValueError("bias must be a float strictly between 0 and 1")

    rng = utils._coerce_rng(rng)

    # Rejection sampling; almost all Boolean functions are non-degenerate
    while True:
        f = random_function_with_bias(n, bias=bias, rng=rng)
        if not f.is_degenerate():
            return f


def random_degenerate_function(
    n: int,
    bias: float = 0.5,
    *,
    rng=None,
) -> BooleanFunction:
    """
    Generate a random degenerate Boolean function.

    A Boolean function is degenerate if at least one variable is
    non-essential, i.e., the function does not depend on that variable.
    This function constructs a degenerate Boolean function by selecting
    one variable uniformly at random and enforcing that the output is
    independent of that variable.

    Parameters
    ----------
    n : int
        Number of Boolean variables.
    bias : float, optional
        Probability that a truth-table entry equals 1 for the underlying
        (n−1)-variable function. Default is 0.5.
    rng : int, np.random.Generator, np.random.RandomState, random.Random, or None, optional
        Random number generator or seed specification. Passed to
        ``utils._coerce_rng``.

    Returns
    -------
    BooleanFunction
        Random degenerate Boolean function on ``n`` variables.

    Raises
    ------
    ValueError
        If ``n`` is not a positive integer.
    ValueError
        If ``bias`` is not strictly between 0 and 1.

    Notes
    -----
    - Exactly one variable is forced to be non-essential by construction,
      though additional variables may also be non-essential by chance.
    - The degenerate variable is chosen uniformly at random.
    - The resulting distribution is not uniform over all degenerate Boolean
      functions.
    - This construction avoids rejection sampling.
    """
    if not isinstance(n, (int, np.integer)) or n <= 0:
        raise ValueError("n must be a positive integer")

    if not isinstance(bias, (float, np.floating)) or not (0.0 < bias < 1.0):
        raise ValueError("bias must be a float strictly between 0 and 1")

    rng = utils._coerce_rng(rng)

    # Generate an (n-1)-variable Boolean function
    f_original = random_function_with_bias(n - 1, bias=bias, rng=rng)

    # Choose the non-essential variable uniformly at random
    index_non_essential_variable = rng.integers(n)

    f = np.zeros(2**n, dtype=np.uint8)

    # Copy the (n-1)-variable function across both values of the non-essential variable
    block = 2 ** index_non_essential_variable
    indices = (np.arange(2**n) // block) % 2 == 1
    f[indices] = f_original.f
    f[~indices] = f_original.f

    return BooleanFunction._from_f_unchecked(f)


def random_non_canalizing_function(
    n: int,
    bias: float = 0.5,
    *,
    rng=None,
) -> BooleanFunction:
    """
    Generate a random non-canalizing Boolean function.

    A Boolean function is canalizing if there exists at least one variable
    and a value of that variable such that fixing it forces the output of
    the function. This function samples Boolean functions from the
    Bernoulli(bias) ensemble until a non-canalizing function is obtained.

    Parameters
    ----------
    n : int
        Number of Boolean variables. Must satisfy ``n > 1``.
    bias : float, optional
        Probability that a truth-table entry equals 1. Default is 0.5.
    rng : int, np.random.Generator, np.random.RandomState, random.Random, or None, optional
        Random number generator or seed specification. Passed to
        ``utils._coerce_rng``.

    Returns
    -------
    BooleanFunction
        Random non-canalizing Boolean function on ``n`` variables.

    Raises
    ------
    ValueError
        If ``n`` is not an integer greater than 1.
    ValueError
        If ``bias`` is not strictly between 0 and 1.

    Notes
    -----
    - This function uses rejection sampling.
    - For moderate bias values, almost all Boolean functions are
      non-canalizing.
    - Extremely biased functions are more likely to be canalizing and may
      lead to longer sampling times.

    References
    ----------
    C. Kadelka, J. Kuipers, and R. Laubenbacher (2017).
    The influence of canalization on the robustness of Boolean networks.
    Physica D: Nonlinear Phenomena, 353, 39–47.
    """
    if not isinstance(n, (int, np.integer)) or n <= 1:
        raise ValueError("n must be an integer greater than 1")

    if not isinstance(bias, (float, np.floating)) or not (0.0 < bias < 1.0):
        raise ValueError("bias must be a float strictly between 0 and 1")

    rng = utils._coerce_rng(rng)

    # Rejection sampling; most Boolean functions are non-canalizing
    while True:
        f = random_function_with_bias(n, bias=bias, rng=rng)
        if not f.is_canalizing():
            return f



def random_non_canalizing_non_degenerate_function(
    n: int,
    bias: float = 0.5,
    *,
    rng=None,
) -> BooleanFunction:
    """
    Generate a random Boolean function that is both non-canalizing and
    non-degenerate.

    A Boolean function is non-canalizing if no variable can force the output
    when fixed, and non-degenerate if every variable is essential. This
    function samples Boolean functions from the Bernoulli(bias) ensemble
    until both properties are satisfied.

    Parameters
    ----------
    n : int
        Number of Boolean variables. Must satisfy ``n > 1``.
    bias : float, optional
        Probability that a truth-table entry equals 1. Default is 0.5.
    rng : int, np.random.Generator, np.random.RandomState, random.Random, or None, optional
        Random number generator or seed specification. Passed to
        ``utils._coerce_rng``.

    Returns
    -------
    BooleanFunction
        Random Boolean function on ``n`` variables that is both non-canalizing
        and non-degenerate.

    Raises
    ------
    ValueError
        If ``n`` is not an integer greater than 1.
    ValueError
        If ``bias`` is not strictly between 0 and 1.

    Notes
    -----
    - This function uses rejection sampling.
    - For moderate bias values and sufficiently large ``n``, almost all
      Boolean functions are both non-canalizing and non-degenerate.
    - Extremely biased functions are more likely to be canalizing or
      degenerate and may lead to longer sampling times.

    References
    ----------
    C. Kadelka, J. Kuipers, and R. Laubenbacher (2017).
    The influence of canalization on the robustness of Boolean networks.
    Physica D: Nonlinear Phenomena, 353, 39–47.
    """
    if not isinstance(n, (int, np.integer)) or n <= 1:
        raise ValueError("n must be an integer greater than 1")

    if not isinstance(bias, (float, np.floating)) or not (0.0 < bias < 1.0):
        raise ValueError("bias must be a float strictly between 0 and 1")

    rng = utils._coerce_rng(rng)

    # Rejection sampling; almost all Boolean functions satisfy both properties
    while True:
        f = random_function_with_bias(n, bias=bias, rng=rng)
        if not f.is_canalizing() and not f.is_degenerate():
            return f

_uniform_over_functions_weights = {}

def _get_uniform_over_functions_weights(max_n, is_ncf=True):
    """
    Compute dynamic-programming weights for uniform canalized layer structures.

    This function constructs a dynamic-programming table ``W`` used to sample
    canalized output bitstrings with probability proportional to the inverse
    factorials of their layer sizes. The table encodes the total weight of all
    valid completions of a partially constructed canalized structure.

    Parameters
    ----------
    max_n : int
        Maximum total length of the canalized output string.
    is_ncf : bool, optional
        If True (default), enforce the nested canalizing function (NCF)
        constraint that the final canalizing layer has size at least 2.
        If False, no constraint is imposed on the final layer size.

    Returns
    -------
    W : ndarray of shape (max_n + 1, max_n + 2)
        Dynamic-programming weight table. ``W[m, s]`` gives the total weight of
        all valid completions with ``m`` positions remaining and current layer
        size ``s``.

    Notes
    -----
    The recursion is given by

        W[m, s] = W[m - 1, s + 1] + (1 / s!) * W[m - 1, 1],

    corresponding to either extending the current canalizing layer or
    terminating it and starting a new layer of size 1. The base case ``m = 0``
    accounts for the weight of the final layer and enforces the optional NCF
    constraint.

    Results are cached internally to avoid recomputation for repeated calls
    with the same ``max_n`` and ``NCF`` values.
    """
    if (max_n,is_ncf) in _uniform_over_functions_weights:
        return _uniform_over_functions_weights[(max_n,is_ncf)]
    
    W = np.zeros((max_n + 1, max_n + 2), dtype=float)

    # Base case: no positions left -> close final run
    inv_factorial_of_s = 1.0
    for s in range(1, max_n + 2):
        inv_factorial_of_s /= s
        if (not is_ncf) or (s >= 2):
            W[0, s] = inv_factorial_of_s
        else:
            W[0, s] = 0.0

    # DP
    for m in range(1, max_n + 1):
        inv_factorial_of_s = 1.0
        for s in range(1, max_n + 1):
            inv_factorial_of_s /= s
            W[m, s] = (
                W[m - 1, s + 1] +
                inv_factorial_of_s * W[m - 1, 1]
            )
    
    _uniform_over_functions_weights[(max_n,is_ncf)] = W    
    return W


def sample_canalized_outputs_uniform_over_functions(n, W, *, rng):
    """
    Sample a canalized output bitstring yielding uniform layer-structure weighting.

    This function samples a binary vector ``b`` of length ``n`` representing
    canalized output values, where consecutive equal values are part of the same
    canalizing layer. The probability of a given layer structure
    ``(k_1, k_2, ..., k_r)`` is proportional to

        1 / (k_1! k_2! ... k_r!).
        
    Sampling is performed sequentially using precomputed dynamic-programming
    weights ``W``, stored in _uniform_over_functions_weights.

    Parameters
    ----------
    n : int
        Length of the output bitstring to sample.
    W : ndarray of shape (n+1, n+1)
        Dynamic-programming weight table, where ``W[m, s]`` gives the total
        weight of all valid completions with ``m`` positions remaining and
        current layer size ``s``.
    rng : numpy.random.Generator
        Random number generator used for sampling.

    Returns
    -------
    b : ndarray of shape (n,), dtype int
        Sampled binary canalized output vector.

    Notes
    -----
    The bitstring is generated left-to-right. At each step, the algorithm
    probabilistically chooses whether to extend the current layer or start a
    new one, using the weights in ``W`` to ensure correct global sampling
    probabilities. The first bit is chosen uniformly at random.
    """
    rng = utils._coerce_rng(rng)

    b = np.zeros(n, dtype=int)

    # Randomly pick first canalized output
    b[0] = rng.integers(2)

    s = 1              # current layer size
    inv_factorial_of_s = 1.0
    m = n - 1          # positions remaining to fill

    for i in range(1, n):

        # DP-weighted decision
        w_extend = W[m - 1, s + 1]
        w_split  = inv_factorial_of_s * W[m - 1, 1]

        p_extend = w_extend / (w_extend + w_split)
        if rng.random() < p_extend:
            b[i] = b[i - 1]
            s += 1
            inv_factorial_of_s /= s
        else:
            b[i] = 1 - b[i - 1]
            s = 1
            inv_factorial_of_s = 1.0

        m -= 1

    return b

def random_k_canalizing_function(
    n: int,
    k: int,
    exact_depth: bool = False,
    uniform_over_functions: bool = True,
    *,
    rng=None,
) -> BooleanFunction:
    """
    Generate a random k-canalizing Boolean function in n variables.

    A Boolean function is k-canalizing if it has at least k conditionally 
    canalizing variables. If ``exact_depth`` is True, the function has exactly
    k conditionally canalizing variables; otherwise, its canalizing depth 
    may exceed k.

    Parameters
    ----------
    n : int
        Number of Boolean variables.
    k : int
        Requested canalizing depth. Must satisfy ``0 <= k <= n``.
        Setting ``k = n`` generates a nested canalizing function.
    exact_depth : bool, optional
        If True, enforce that the canalizing depth is exactly ``k``.
        If False (default), the depth is at least ``k``.
    uniform_over_functions : bool, optional
        If True (default), the function is sampled uniformly at random
        from the set of Boolean functions consistent with the specified
        constraints (n, k, exact_depth).
        
        Internally, this is achieved by a rejection sampling scheme that
        compensates for all combinatorial multiplicities arising from
        symmetric canalizing layers and possible merges between the
        outer canalizing layers and the core function.
    
        If False, canalized outputs are sampled independently as a bitstring,
        and no rejection correction is applied. In this case, the resulting
        distribution is biased toward Boolean functions with more symmetric
        canalizing structures.
    rng : int, numpy.random.Generator, numpy.random.RandomState, random.Random, or None, optional
        Random number generator or seed specification. Passed to
        ``utils._coerce_rng``.

    Returns
    -------
    BooleanFunction
        A Boolean function on ``n`` variables with canalizing depth at
        least ``k`` (or exactly ``k`` if ``exact_depth=True``).

    Raises
    ------
    AssertionError
        If ``n`` is not a positive integer.
    AssertionError
        If ``k`` does not satisfy ``0 <= k <= n``.
    AssertionError
        If ``exact_depth=True`` and ``k = n-1`` (no such functions exist).

    Notes
    -----
    When uniform_over_functions=True, this function samples uniformly from the
    space of Boolean functions with the specified canalizing properties.
    As a consequence, different canalizing layer structures generally
    occur with different frequencies, reflecting the fact that they
    support different numbers of Boolean functions.
    
    Uniformity over canalizing layer structures is not enforced and is
    not expected.

    The construction follows the standard decomposition of a k-canalizing
    function into canalizing variables, canalizing inputs and outputs, and
    a residual core function on ``n-k`` variables.

    References
    ----------
    He, Q., and Macauley, M. (2016).
        Stratification and enumeration of Boolean functions by canalizing depth.
        Physica D: Nonlinear Phenomena, 314, 1–8.

    Dimitrova, E., Stigler, B., Kadelka, C., and Murrugarra, D. (2022).
        Revealing the canalizing structure of Boolean functions: Algorithms
        and applications. Automatica, 146, 110630.
    """
    rng = utils._coerce_rng(rng)

    assert isinstance(n, (int, np.integer)) and n > 0, "n must be a positive integer"
    assert n - k != 1 or not exact_depth, (
        "There are no functions of exact canalizing depth n-1.\nEither set exact_depth=False or ensure k != n-1"
    )
    assert isinstance(k, (int, np.integer)) and 0 <= k and k <= n, (
        "k, the canalizing depth, must satisfy 0 <= k <= n."
    )
    
    if k==n-1 and n>1: #canalizing functions with depth n-1>0 really have depth n 
        k=n

    # Step 1: canalizing inputs and variables
    aas = rng.integers(2, size=k)
    can_vars = rng.choice(n, k, replace=False)

    while True: 
        #include the generation of canalized outputs in the rejection sampling scheme
        #because some output vectors `bbs` may give rise to more distinct functions than others
        
        # Step 2: canalized outputs, determining layers
        if uniform_over_functions:
            bbs = sample_canalized_outputs_uniform_over_functions(
                k,
                _get_uniform_over_functions_weights(k, is_ncf=(k >= n)),
                rng=rng,
            )
        else:
            bbs = rng.integers(2, size=k)
            
        # Step 3: sample core function using efficient rejection sampling
        if k < n:
            core_function = random_function(
                n=n - k,
                depth=0,
                exact_depth=exact_depth,#True if exact_depth or uniform_over_functions else False,
                allow_degenerate_functions=False,
                rng=rng,
            )
            
            if exact_depth or k==0 or not uniform_over_functions:
                break
            else:
                #check if the core function is canalizing and correct for combinatorial
                #bias in the selection of the final functions
                
                #compute canalizing info of core_function, stored in f.properties
                core_function.get_layer_structure() 
                if core_function.properties['CanalizingDepth'] == 0:
                    break
                else: #rejection sampling is efficient because it happens with probability <= 50%
                    bbs_core = core_function.properties['CanalizedOutputs']
                    if bbs_core[0] != bbs[-1]:  
                        #no combinatorial explosion, the core varaibles start a new layer
                        break
                    else: 
                        # merge occurs: last layer grows
                        
                        # determine size of last outer layer in bbs
                        s = 1
                        for i in range(k - 1, 0, -1):
                            if bbs[i] == bbs[i - 1]:
                                s += 1
                            else:
                                break
                            
                        # determine number of additional variables in the same layer
                        s_core = 1
                        for i in range(core_function.properties['CanalizingDepth'] - 1):
                            if bbs_core[i] == bbs_core[i + 1]:
                                s_core += 1
                            else:
                                break                    
                        
                        accept_prob = 1 / math.comb(s + s_core, s)
                        
                        if rng.random() <= accept_prob:
                            break
        else:
            core_function = [1 - bbs[-1]]
            break
        
    # Step 4: build truth table and return canalizing Boolean function
    left_side_of_truth_table = utils.get_left_side_of_truth_table(n)
    f = np.full(2**n, -1, dtype=np.int8)
    
    for j in range(k):
        mask = (left_side_of_truth_table[:, can_vars[j]] == aas[j]) & (f < 0)
        f[mask] = bbs[j]
        
    # fill remaining with core truth table
    f[f < 0] = np.asarray(core_function, dtype=np.int8)

    return BooleanFunction._from_f_unchecked(f)


def random_k_canalizing_function_with_specific_layer_structure(
    n: int,
    layer_structure: list,
    exact_depth: bool = False,
    *,
    rng=None,
) -> BooleanFunction:
    """
    Generate a random Boolean function with a specified canalizing layer structure.

    The canalizing layer structure is given as a list
    ``[k_1, ..., k_r]``, where each ``k_i`` specifies the number of
    canalizing variables in the i-th layer. The total canalizing depth is
    ``sum(layer_structure)``.

    If ``sum(layer_structure) == n`` and ``n > 1``, the function is a nested
    canalizing function and the final layer is required to have size at
    least 2.

    Parameters
    ----------
    n : int
        Total number of Boolean variables.
    layer_structure : list of int
        Canalizing layer structure ``[k_1, ..., k_r]``. Each entry must be
        at least 1. If ``sum(layer_structure) == n`` and ``n > 1``, the final
        entry must satisfy ``layer_structure[-1] >= 2``.
    exact_depth : bool, optional
        If True, enforce that the canalizing depth is exactly
        ``sum(layer_structure)``. If False (default), additional canalizing
        variables may occur in the core function.
    rng : int, numpy.random.Generator, numpy.random.RandomState, random.Random, or None, optional
        Random number generator or seed specification. Passed to
        ``utils._coerce_rng``.

    Returns
    -------
    BooleanFunction
        A Boolean function on ``n`` variables with the prescribed canalizing
        layer structure, plus potentially more canalizing variables (if exact_depth=True).

    Raises
    ------
    AssertionError
        If ``n`` is not a positive integer.
    AssertionError
        If ``sum(layer_structure)`` does not satisfy ``0 <= sum(layer_structure) <= n``.
    AssertionError
        If ``exact_depth=True`` and ``sum(layer_structure) = n - 1``.
    AssertionError
        If ``sum(layer_structure) = n > 1`` and the final layer has size less
        than 2.
    AssertionError
        If any entry of ``layer_structure`` is less than 1.

    Notes
    -----
    For fixed parameter values, this function samples uniformly at random
    from the ensemble of Boolean functions consistent with the specified
    canalizing layer structure and additional constraints.

    The construction follows the standard decomposition of a canalizing
    function into ordered canalizing layers and a residual core function on
    the remaining variables.

    References
    ----------
    He, Q., and Macauley, M. (2016).
        Stratification and enumeration of Boolean functions by canalizing depth.
        Physica D: Nonlinear Phenomena, 314, 1–8.

    Kadelka, C., Kuipers, J., and Laubenbacher, R. (2017).
        The influence of canalization on the robustness of Boolean networks.
        Physica D: Nonlinear Phenomena, 353, 39–47.
    """
    rng = utils._coerce_rng(rng)
    depth = sum(layer_structure)  # canalizing depth
    if depth == 0:
        layer_structure = [0]

    assert isinstance(n, (int, np.integer)) and n > 0, "n must be an integer > 0"
    assert n - depth != 1 or not exact_depth, (
        "There are no functions of exact canalizing depth n-1.\nEither set exact_depth=False or ensure depth=sum(layer_structure)!=n-1."
    )
    assert 0 <= depth and depth <= n, "Ensure 0 <= depth = sum(layer_structure) <= n."
    assert depth < n or layer_structure[-1] > 1 or n == 1, (
        "The last layer of an NCF (i.e., an n-canalizing function) has to have size >= 2 whenever n > 1.\nIf depth=sum(layer_structure)=n, ensure that layer_structure[-1]>=2."
    )
    assert min(layer_structure) >= 1, (
        "Each layer must have at least one variable (each element of layer_structure must be >= 1)."
    )

    size_state_space = 2**n
    aas = rng.integers(2, size=depth)  # canalizing inputs
    b0 = rng.integers(2)
    bbs = [b0] * layer_structure[0]  # canalized outputs for first layer
    for i in range(1, len(layer_structure)):
        if i % 2 == 0:
            bbs.extend([b0] * layer_structure[i])
        else:
            bbs.extend([1 - b0] * layer_structure[i])
    can_vars = rng.choice(n, depth, replace=False)
    f = np.zeros(size_state_space, dtype=int)
    if depth < n:
        core_function = random_function(
            n=n - depth,
            depth=0,
            exact_depth=exact_depth,
            allow_degenerate_functions=False,
            rng=rng,
        )
    else:
        core_function = [1 - bbs[-1]]

    left_side_of_truth_table = utils.get_left_side_of_truth_table(n)
    f = np.full(2**n, -1, dtype=np.int8)
    for j in range(depth):
        mask = (left_side_of_truth_table[:, can_vars[j]] == aas[j]) & (f < 0)
        f[mask] = bbs[j]
    # fill remaining with core truth table
    f[f < 0] = np.asarray(core_function, dtype=np.int8)

    return BooleanFunction._from_f_unchecked(f)


def random_NCF(
    n: int,
    uniform_over_functions: bool = True,
    layer_structure: list | None = None,
    *,
    rng=None
) -> BooleanFunction:
    """
    Generate a random nested canalizing Boolean function in n variables.

    A nested canalizing function (NCF) is an n-canalizing Boolean function,
    i.e., a function whose canalizing depth equals the number of variables.
    Optionally, a specific canalizing layer structure may be prescribed.

    Parameters
    ----------
    n : int
        Total number of Boolean variables.
    uniform_over_functions : bool, optional
        If True (default) and ``layer_structure`` is None, canalizing layer
        structures are sampled uniformly at random, removing the bias toward
        symmetric structures induced by independent sampling of canalized
        outputs. If False, canalized outputs are sampled independently and
        uniformly as a bitstring, which biases the distribution toward more
        symmetric layer structures.
        This parameter is ignored if ``layer_structure`` is provided.
    layer_structure : list of int or None, optional
        Canalizing layer structure ``[k_1, ..., k_r]``. Each entry must be at
        least 1. If provided, it must satisfy ``sum(layer_structure) == n``.
        If ``n > 1``, the final entry must satisfy
        ``layer_structure[-1] >= 2``. If None (default), the layer structure
        is sampled at random.
    rng : int, numpy.random.Generator, numpy.random.RandomState, random.Random, or None, optional
        Random number generator or seed specification. Passed to
        ``utils._coerce_rng``.

    Returns
    -------
    BooleanFunction
        A nested canalizing Boolean function on ``n`` variables.

    Raises
    ------
    AssertionError
        If ``n`` is not a positive integer.
    AssertionError
        If ``layer_structure`` is provided but does not satisfy
        ``sum(layer_structure) == n``.
    AssertionError
        If ``n > 1`` and the final layer has size less than 2.

    Notes
    -----
    For fixed parameter values, this function samples uniformly at random
    from the ensemble of nested canalizing Boolean functions consistent with
    the specified constraints. Non-uniformity arises only when
    ``uniform_over_functions=False``.

    This function is a convenience wrapper around
    ``random_k_canalizing_function`` and
    ``random_k_canalizing_function_with_specific_layer_structure``.

    References
    ----------
    He, Q., and Macauley, M. (2016).
        Stratification and enumeration of Boolean functions by canalizing depth.
        Physica D: Nonlinear Phenomena, 314, 1–8.

    Kadelka, C., Kuipers, J., and Laubenbacher, R. (2017).
        The influence of canalization on the robustness of Boolean networks.
        Physica D: Nonlinear Phenomena, 353, 39–47.
    """
    rng = utils._coerce_rng(rng)
    if layer_structure is None:
        return random_k_canalizing_function(n, 
                                            n, 
                                            uniform_over_functions=uniform_over_functions,
                                            exact_depth=False, 
                                            rng=rng)
    else:
        assert sum(layer_structure) == n, "Ensure sum(layer_structure) == n."
        assert layer_structure[-1] > 1 or n == 1, (
            "The last layer of an NCF has to have size >= 2 whenever n > 1.\nEnsure that layer_structure[-1]>=2."
        )
        return random_k_canalizing_function_with_specific_layer_structure(
            n, 
            layer_structure, 
            exact_depth=False, 
            rng=rng
        )


## Random network generation
def random_degrees(
    N: int,
    n: int | float | list | np.ndarray,
    indegree_distribution: str = "constant",
    allow_self_loops: bool = False,
    *,
    rng=None,
) -> np.ndarray:
    """
    Draw an in-degree vector for a directed network with N nodes.

    This function either accepts a user-specified in-degree vector or
    samples in-degrees independently for each node from a specified
    distribution.

    Parameters
    ----------
    N : int
        Number of nodes in the network. Must be a positive integer.
    n : int, float, list of int, or ndarray of int
        Interpretation depends on ``indegree_distribution``:

        - If ``n`` is a length-``N`` vector of integers, it is interpreted
          as a user-specified in-degree sequence and returned after
          validation.

        - If ``indegree_distribution`` is one of
          ``{'constant', 'dirac', 'delta'}``, then ``n`` is a single integer
          specifying the in-degree of every node.

        - If ``indegree_distribution == 'uniform'``, then ``n`` is a positive
          integer upper bound, and each node independently receives an
          in-degree sampled uniformly from
          ``{1, 2, ..., n}``.

        - If ``indegree_distribution == 'poisson'``, then ``n`` is the Poisson
          rate parameter ``λ > 0``. Each node independently receives a
          Poisson(``λ``) draw, truncated to lie in
          ``[1, N - int(not allow_self_loops)]``.

    indegree_distribution : str, optional
        Distribution used to generate in-degrees when ``n`` is not a vector.
        Must be one of ``{'constant', 'dirac', 'delta', 'uniform', 'poisson'}``.
        Default is ``'constant'``.
    allow_self_loops : bool, optional
        If True, in-degrees may be as large as ``N``.
        If False (default), self-loops are disallowed in subsequent wiring
        generation. This is enforced here by capping in-degrees at ``N-1``.
    rng : int, numpy.random.Generator, numpy.random.RandomState, random.Random, or None, optional
        Random number generator or seed specification. Passed to
        ``utils._coerce_rng``.

    Returns
    -------
    indegrees : ndarray of int, shape (N,)
        In-degree of each node. For sampled distributions, values lie in
        ``[1, N - int(not allow_self_loops)]``.

    Raises
    ------
    AssertionError
        If inputs are malformed, out of range, or an unsupported distribution
        is requested.

    Notes
    -----
    When sampling is requested, in-degrees for different nodes are generated
    independently. No attempt is made to enforce graphicality or feasibility
    of the resulting degree sequence for a particular wiring model; such
    constraints must be handled downstream.

    Examples
    --------
    >>> random_degrees(5, n=2, indegree_distribution='constant')
    array([2, 2, 2, 2, 2])

    >>> random_degrees(4, n=2, indegree_distribution='uniform', allow_self_loops=False)
    array([2, 1, 2, 2])

    >>> random_degrees(6, n=1.7, indegree_distribution='poisson')
    array([1, 2, 1, 1, 2, 1])

    >>> random_degrees(3, n=[1, 2, 1])
    array([1, 2, 1])
    """
    rng = utils._coerce_rng(rng)

    if isinstance(n, (list, np.ndarray)):
        assert (
            utils.is_list_or_array_of_ints(n, required_length=N)
            and min(n) >= 1
            and max(n) <= N - int(not allow_self_loops)
        ), (
            "A vector n was submitted.\nEnsure that n is an N-dimensional vector where each element is an integer between 1 and "
            + ("N-1" if not allow_self_loops else "N")
            + " representing the indegree of each nodde."
        )
        indegrees = np.array(n, dtype=int)
    elif indegree_distribution.lower() in ["constant", "dirac", "delta"]:
        assert (
            isinstance(n, (int, np.integer))
            and n >= 1
            and n <= N - int(not allow_self_loops)
        ), (
            "n must be an integer between 1 and "
            + ("N-1" if not allow_self_loops else "N")
            + " describing the constant degree of each node."
        )
        indegrees = np.ones(N, dtype=int) * n
    elif indegree_distribution.lower() == "uniform":
        assert (
            isinstance(n, (int, np.integer))
            and n >= 1
            and n <= N - int(not allow_self_loops)
        ), (
            "n must be an integer between 1 and "
            + ("N-1" if not allow_self_loops else "N")
            + " representing the upper bound of a uniform degree distribution (lower bound == 1)."
        )
        indegrees = rng.integers(1, n + 1, size=N)
    elif indegree_distribution.lower() == "poisson":
        assert isinstance(n, (int, float, np.integer, np.floating)) and n > 0, (
            "n must be a float > 0 representing the Poisson parameter."
        )
        indegrees = np.maximum(
            np.minimum(rng.poisson(lam=n, size=N), N - int(not allow_self_loops)), 1
        )
    else:
        raise AssertionError(
            "None of the predefined in-degree distributions were chosen.\nTo use a user-defined in-degree vector, submit an N-dimensional vector as argument for n; each element of n must an integer between 1 and N."
        )
    return indegrees


def random_edge_list(
    N: int,
    indegrees: Sequence[int],
    allow_self_loops: bool,
    min_out_degree_one: bool = False,
    *,
    rng=None,
) -> list:
    """
    Generate a random directed edge list for a network with prescribed in-degrees.

    Each node ``i`` receives exactly ``indegrees[i]`` incoming edges, with
    regulators chosen uniformly at random from the set of admissible source
    nodes. Optionally, the construction enforces that every node regulates at
    least one other node.

    Parameters
    ----------
    N : int
        Number of nodes in the network.
    indegrees : sequence of int
        Length-``N`` sequence specifying the number of incoming edges for each
        node.
    allow_self_loops : bool
        If True, self-loops (edges from a node to itself) are allowed.
        Default is False.
    min_out_degree_one : bool, optional
        If True, enforce that every node has at least one outgoing edge.
        This is achieved by rewiring edges while preserving the prescribed
        in-degree sequence. Default is False.
    rng : int, numpy.random.Generator, numpy.random.RandomState, random.Random, or None, optional
        Random number generator or seed specification. Passed to
        ``utils._coerce_rng``.

    Returns
    -------
    edge_list : list of tuple of int
        List of directed edges represented as ``(source, target)`` pairs.

    Raises
    ------
    ValueError
        If ``N`` or ``indegrees`` are inconsistent.
    AssertionError
        If sampling constraints cannot be satisfied.

    Notes
    -----
    Regulators for each node are sampled uniformly at random without
    replacement from the set of admissible source nodes. If
    ``min_out_degree_one``, the algorithm post-processes
    the initially sampled edge list by replacing edges until every node
    has at least one outgoing edge, while preserving all in-degrees and
    respecting the self-regulation constraint.

    No guarantee is made that the resulting edge list is uniformly sampled
    from the space of all directed graphs satisfying the constraints.
    """

    rng = utils._coerce_rng(rng)

    # ------------------------------------------------------------
    # Step 1: generate initial edge list
    # ------------------------------------------------------------
    edge_list = []
    for i in range(N):
        if not allow_self_loops:
            candidates = np.append(np.arange(i), np.arange(i + 1, N))
        else:
            candidates = np.arange(N)

        indices = rng.choice(candidates, indegrees[i], replace=False)
        edge_list.extend(zip(indices, np.full(indegrees[i], i, dtype=int)))

    # ------------------------------------------------------------
    # Step 2: enforce at least one outgoing edge per node (optional)
    # ------------------------------------------------------------
    if min_out_degree_one:
        target_sources = [set() for _ in range(N)]
        outdegrees = np.zeros(N, dtype=int)

        for s, t in edge_list:
            target_sources[t].add(s)
            outdegrees[s] += 1

        sum_indegrees = len(edge_list)

        while np.min(outdegrees) == 0:
            index_sink = np.where(outdegrees == 0)[0][0]
            index_edge = rng.integers(sum_indegrees)

            old_source, t = edge_list[index_edge]

            if not allow_self_loops and t == index_sink:
                continue
            if index_sink in target_sources[t]:
                continue

            # perform replacement
            target_sources[t].discard(old_source)
            target_sources[t].add(index_sink)

            edge_list[index_edge] = (index_sink, t)

            outdegrees[index_sink] += 1
            outdegrees[old_source] -= 1

    return edge_list


def random_wiring_diagram(
    N: int,
    n: int | float | list | np.ndarray,
    allow_self_loops: bool = False,
    strongly_connected: bool = False,
    indegree_distribution: str = "constant",
    min_out_degree_one: bool = False,
    max_strong_connectivity_attempts: int = 1000,
    *,
    rng=None,
) -> tuple:
    """
    Generate a random wiring diagram for a directed network with N nodes.

    A wiring diagram specifies, for each node, the set of its regulators
    (incoming neighbors). In-degrees are first generated according to the
    specified distribution, after which edges are sampled uniformly at
    random subject to the requested constraints.

    Parameters
    ----------
    N : int
        Number of nodes in the network.
    n : int, float, list of int, or ndarray of int
        Parameter determining the in-degree sequence. Interpretation depends
        on ``indegree_distribution``:

        - If a length-``N`` vector is provided, it is interpreted as the
          in-degree of each node.
        - If ``indegree_distribution`` is ``'constant'`` (or ``'dirac'`` /
          ``'delta'``), ``n`` specifies the in-degree of every node.
        - If ``indegree_distribution`` is ``'uniform'``, ``n`` specifies the
          upper bound of a discrete uniform distribution on
          ``{1, ..., n}``.
        - If ``indegree_distribution`` is ``'poisson'``, ``n`` is the Poisson
          rate parameter ``lambda > 0``.
    allow_self_loops : bool, optional
        If True, self-loops (edges from a node to itself) are allowed.
        Default is False.
    strongly_connected : bool, optional
        If True, repeatedly resample the wiring diagram until a strongly
        connected network is obtained, or until the maximum number of
        attempts is exceeded. Default is False.
    indegree_distribution : str, optional
        Distribution used to generate in-degrees. Must be one of
        ``{'constant', 'dirac', 'delta', 'uniform', 'poisson'}``.
        Default is ``'constant'``.
    min_out_degree_one : bool, optional
        If True, enforce that every node has at least one outgoing edge.
        This is achieved by rewiring edges while preserving the in-degree
        sequence. Default is False.
    max_strong_connectivity_attempts : int, optional
        Maximum number of attempts to generate a strongly connected wiring
        diagram before raising a ``RuntimeError``. Default is 1000.
    rng : int, numpy.random.Generator, numpy.random.RandomState, random.Random, or None, optional
        Random number generator or seed specification. Passed to
        ``utils._coerce_rng``.

    Returns
    -------
    WiringDiagram
        A wiring diagram object encoding the regulator set of each node.

    Raises
    ------
    RuntimeError
        If ``strongly_connected=True`` and a strongly connected wiring diagram
        cannot be generated within the specified number of attempts.

    Notes
    -----
    In-degrees are generated first using ``random_degrees``, and edges are
    then sampled uniformly at random subject to the imposed constraints.
    When ``strongly_connected=True`` or
    ``min_out_degree_one=True``, the resulting distribution is
    not uniform over all wiring diagrams with the given in-degree sequence.

    This function is a high-level convenience wrapper around
    ``random_degrees`` and ``random_edge_list``.

    Examples
    --------
    >>> W = random_wiring_diagram(5, n=2)
    >>> W = random_wiring_diagram(10, n=3, strongly_connected=True)
    >>> W = random_wiring_diagram(6, n=[1, 2, 1, 2, 1, 2])
    """
    rng = utils._coerce_rng(rng)
    indegrees = random_degrees(
        N,
        n,
        indegree_distribution=indegree_distribution,
        allow_self_loops=allow_self_loops,
        rng=rng,
    )

    counter = 0
    while True:  # Keep generating until we have a strongly connected graph
        edges_wiring_diagram = random_edge_list(
            N,
            indegrees,
            allow_self_loops,
            min_out_degree_one=min_out_degree_one,
            rng=rng,
        )
        if strongly_connected:
            # may take a long time ("forever") if n is small and N is large
            G = nx.from_edgelist(edges_wiring_diagram, create_using=nx.MultiDiGraph())
            if not nx.is_strongly_connected(G):
                counter += 1
                if counter > max_strong_connectivity_attempts:
                    raise RuntimeError(
                        "Made "
                        + str(max_strong_connectivity_attempts)
                        + " unsuccessful attempts to generate a strongly connected wiring diagram of "
                        + str(N)
                        + " nodes and degrees "
                        + str(indegrees)
                        + ".\nYou may increase the number of attempts by modulating the parameter max_strong_connectivity_attempts."
                    )
                continue
        break
    I = [[] for _ in range(N)]
    for edge in edges_wiring_diagram:
        I[edge[1]].append(edge[0])
    for i in range(N):
        I[i] = np.sort(I[i])
    return WiringDiagram(I)


def rewire_wiring_diagram(
    I: list | np.ndarray | WiringDiagram,
    average_swaps_per_edge: float = 50,
    allow_new_self_loops: bool = False,
    allow_self_loop_rewiring: bool = False,
    *,
    rng=None,
) -> list:
    """
    Degree-preserving rewiring of a wiring diagram via double-edge swaps.

    The wiring diagram is represented in regulator form: ``I[target]`` lists
    all regulators (incoming neighbors) of ``target``. The algorithm performs
    random double-edge swaps of the form
    ``(u -> v, x -> y) -> (u -> y, x -> v)``, while preserving both the
    in-degree and out-degree of every node. Parallel edges are disallowed.

    Parameters
    ----------
    I : list of array-like or WiringDiagram
        Wiring diagram in regulator representation. For each node ``target``,
        ``I[target]`` contains the regulators of that node. Regulator indices
        must be integers in ``{0, ..., N-1}``. If a ``WiringDiagram`` is
        provided, its internal adjacency representation is used.
    average_swaps_per_edge : float, optional
        Target number of successful double-edge swaps per edge. Larger values
        typically yield better mixing but increase runtime. Default is 50.
    allow_new_self_loops : bool, optional
        If True, new self-loops may be introduced. If False (default), proposed 
        swaps that would introduce a new self-loop are rejected. 
    allow_self_loop_rewiring : bool, optional
        If True, existing self-loops may be rewired. If False (default), 
        existing self-loops are kept fixed and excluded from the pool of 
        swappable edges. 
    rng : int, numpy.random.Generator, numpy.random.RandomState, random.Random, or None, optional
        Random number generator or seed specification. Passed to
        ``utils._coerce_rng``.

    Returns
    -------
    WiringDiagram
        A new wiring diagram obtained by degree-preserving rewiring of ``I``.

    Raises
    ------
    ValueError
        If the input wiring diagram is malformed.
    AssertionError
        If rewiring constraints cannot be satisfied.

    Notes
    -----
    Both in-degrees and out-degrees of all nodes are preserved exactly.
    Duplicate edges are never introduced. Control over self-regulation is
    governed by the two Boolean flags above.

    The resulting wiring diagram is not guaranteed to be sampled uniformly
    from the space of all directed graphs with the same degree sequence; the
    procedure is intended as a practical degree-preserving randomization
    method rather than an exact uniform sampler.

    Examples
    --------
    >>> I = random_network(8,3)
    >>> J = rewire_wiring_diagram(I)
    >>> I.indegrees == J.indegrees
    True
    >>> I.get_outdegrees() == J.get_outdegrees()
    True
    """
    rng = utils._coerce_rng(rng)
    
    if isinstance(I, WiringDiagram):
        I = I.I
    
    N = len(I)

    edges = [
        (int(regulator), target)
        for target in range(N)
        for regulator in I[target]
        if regulator != target or allow_self_loop_rewiring
    ]
    n_total_edges = len(edges)

    Jset = [set(regulators) for regulators in I]

    n_rewires_before_stop = int(average_swaps_per_edge * n_total_edges)
    successes = 0
    attempts = 0
    max_attempts = 50 * n_rewires_before_stop + 100

    # Helper to check if adding edge (regulator->target) is allowed
    def edge_ok(regulator, target):
        if not allow_new_self_loops and regulator == target:
            return False
        if regulator in Jset[target]:
            return False
        return True

    while successes < n_rewires_before_stop and attempts < max_attempts:
        attempts += 1

        # Pick two distinct edges uniformly at random
        i, j = rng.choice(n_total_edges, 2, replace=False)

        (u, v) = edges[i]
        (x, y) = edges[j]

        # Swapping identical sources or identical targets is fine in principle,
        # but skip trivial cases that do nothing or re-create the same edges.
        if (u == x) or (v == y):
            continue

        # Proposed swapped edges
        a, b = u, y
        c, d = x, v

        # If the proposed edges are identical to originals, skip
        if (a, b) == (u, v) or (c, d) == (x, y):
            continue

        # Check constraints for both new edges
        if not edge_ok(a, b) or not edge_ok(c, d):
            continue

        # Perform the swap: update adjacency and edge list
        # Remove old edges
        Jset[v].discard(u)
        Jset[y].discard(x)
        # Add new edges
        Jset[b].add(a)
        Jset[d].add(c)
        # Commit edges
        edges[i] = (a, b)
        edges[j] = (c, d)

        successes += 1

    # Reconstruct J from adjacency sets
    J = [np.sort(list(Jset[target])) for target in range(N)]
    return WiringDiagram(J)

def random_network(
    N: int | None = None,
    n: int | float | list | np.ndarray | None = None,
    depth: int | list | np.ndarray = 0,
    exact_depth: bool = False,
    uniform_over_functions: bool = True,
    layer_structure: list | None = None,
    allow_degenerate_functions: bool = False,
    parity: bool = False,
    bias: float | list | np.ndarray = 0.5,
    absolute_bias: float | list | np.ndarray = 0.0,
    use_absolute_bias: bool = True,
    hamming_weight: int | list | np.ndarray | None = None,
    allow_self_loops: bool = False,
    strongly_connected: bool = False,
    indegree_distribution: str = "constant",
    min_out_degree_one: bool = False,
    max_strong_connectivity_attempts: int = 1000,
    I: list | np.ndarray | None | WiringDiagram | nx.DiGraph = None,
    *,
    rng=None,
) -> BooleanNetwork:
    """
    Construct a random Boolean network with configurable wiring and update rules.
    
    The network is built in two stages:
    
    1. Wiring diagram
       If ``I`` is provided, it is used directly as the wiring diagram, where
       ``I[v]`` lists the regulators of node ``v``. Otherwise, a wiring diagram
       for ``N`` nodes is sampled using ``random_wiring_diagram``, with in-degrees
       determined by ``n`` and ``indegree_distribution``. Self-loops may be
       disallowed and strong connectivity may be enforced.
    
    2. Update rules
       For each node ``i``, a Boolean update function with arity
       ``indegrees[i]`` is generated using ``random_function`` subject to the
       requested constraints on canalizing depth or layer structure, parity,
       bias or absolute bias, and exact Hamming weight.
    
    Parameters
    ----------
    N : int or None, optional
        Number of nodes. Required when ``I`` is not provided. Ignored if ``I`` is
        given.
    n : int, float, list of int, ndarray of int, or None, optional
        Controls the in-degree distribution when generating a wiring diagram
        (ignored if ``I`` is given). Interpretation depends on
        ``indegree_distribution``:
    
        - ``'constant'``, ``'dirac'``, ``'delta'``:
          Every node has constant in-degree ``n``.
        - ``'uniform'``:
          ``n`` is an integer upper bound; each node’s in-degree is sampled
          uniformly from ``{1, ..., n}``.
        - ``'poisson'``:
          ``n`` is a positive rate parameter lambda; in-degrees are Poisson(lambda) 
          draws truncated to ``[1, N - int(not allow_self_loops)]``.
        - If ``n`` is a length-``N`` vector of integers, it is taken as the exact
          in-degree sequence.
    depth : int, list of int, or ndarray of int, optional
        Requested canalizing depth per node. If an integer, it is broadcast to
        all nodes and clipped at each node’s in-degree. If a vector, it must have
        length ``N``. Interpreted as a minimum depth unless ``exact_depth=True``.
        Default is 0.
    exact_depth : bool, optional
        If True, each Boolean function is generated with exactly the requested
        canalizing depth (or exactly ``sum(layer_structure[i])`` if a layer
        structure is provided). If False, the canalizing depth is at least as
        large as requested. Default is False.
    uniform_over_functions : bool, optional
        Controls how canalized outputs are sampled when generating canalizing
        functions.
    
        If True (default), canalizing layer structures are sampled uniformly at
        random, i.e., proportional to the inverse factorials of layer sizes,
        removing the bias toward symmetric structures induced by independent
        sampling of canalized outputs.
    
        If False, canalized outputs are sampled independently and uniformly as
        bitstrings, which biases the distribution toward more symmetric layer
        structures.
    
        This parameter is ignored when ``layer_structure`` is explicitly provided.
    layer_structure : list, list of lists, or None, optional
        Canalizing layer structure specifications.
    
        - If None (default), rule generation is controlled by ``depth`` and
          ``exact_depth``.
        - If a single list ``[k1, ..., kr]``, the same structure is used for all
          nodes.
        - If a list of lists of length ``N``, ``layer_structure[i]`` is used for
          node ``i``.
    
        In all cases, ``sum(layer_structure[i])`` must not exceed the in-degree
        of node ``i``. When provided, ``layer_structure`` takes precedence over
        ``depth``.
    allow_degenerate_functions : bool, optional
        If True and ``depth == 0`` and ``layer_structure is None``, degenerate
        Boolean functions (with non-essential inputs) may be generated, as in
        classical NK-Kauffman models. If False, generated functions are required
        to be non-degenerate whenever possible. Default is False.
    parity : bool, optional
        If True, parity Boolean functions are generated for all nodes and all
        other rule parameters are ignored. Default is False.
    bias : float, list of float, or ndarray of float, optional
        Probability of output 1 when generating random (non-canalizing) Boolean
        functions. Used only when ``depth == 0``, ``layer_structure is None``,
        ``parity`` is False, and ``use_absolute_bias`` is False. Scalars are
        broadcast to length ``N``. Must lie in ``[0, 1]``. Default is 0.5.
    absolute_bias : float, list of float, or ndarray of float, optional
        Absolute deviation from 0.5 used when ``use_absolute_bias`` is True.
        Scalars are broadcast to length ``N``. Must lie in ``[0, 1]``. Default 0.0.
    use_absolute_bias : bool, optional
        If True, the bias of each rule is chosen at random from
        ``{0.5*(1-absolute_bias), 0.5*(1+absolute_bias)}``. If False, ``bias`` is
        used directly. Default is True.
    hamming_weight : int, list of int, ndarray of int, or None, optional
        Exact Hamming weight (number of ones) of each truth table. Scalars are
        broadcast to length ``N``. Values must lie in ``{0, ..., 2^k}`` for a
        k-input function. Additional restrictions apply when requesting exact
        depth zero. Default is None.
    allow_self_loops : bool, optional
        If True, self-loops (edges from a node to itself) are allowed.
        Default is False. Ignored if ``I`` is provided.
    strongly_connected : bool, optional
        If True, wiring generation is repeated until a strongly connected
        directed graph is obtained or the attempt limit is exceeded. Ignored if
        ``I`` is provided. Default is False.
    indegree_distribution : str, optional
        Distribution used when sampling in-degrees. Must be one of
        ``{'constant', 'dirac', 'delta', 'uniform', 'poisson'}``. Default
        ``'constant'``.
    min_out_degree_one : bool, optional
        If True, ensure that each node has at least one outgoing edge in the
        generated wiring diagram. Default is False.
    max_strong_connectivity_attempts : int, optional
        Maximum number of attempts to generate a strongly connected wiring
        diagram before raising an error. Default is 1000.
    I : list, ndarray, WiringDiagram, networkx.DiGraph, or None, optional
        Existing wiring diagram. If provided, ``N`` and ``n`` are ignored and
        in-degrees are inferred from ``I``. If I is a BooleanNetwork, its wiring
        diagram is reused and its Boolean update rules are ignored.
    rng : int, numpy.random.Generator, numpy.random.RandomState, random.Random, or None, optional
        Random number generator or seed specification. Passed to
        ``utils._coerce_rng``.
    
    Returns
    -------
    BooleanNetwork
        A Boolean network with wiring diagram ``I`` (given or generated) and
        Boolean update functions generated according to the specified constraints.
    
    Raises
    ------
    AssertionError
        If input shapes or parameter combinations are invalid.
    RuntimeError
        If ``strongly_connected=True`` and a strongly connected wiring diagram
        cannot be generated within the specified number of attempts.
    
    Notes
    -----
    Constraint precedence for rule generation is:
    ``parity`` -> ``layer_structure`` -> ``depth`` / ``exact_depth`` -> bias or
    Hamming-weight constraints.
    
    When ``exact_depth=True`` and the requested depth is zero, Hamming weights
    ``{0, 1, 2^k - 1, 2^k}`` correspond to canalizing functions and are therefore
    disallowed.
    
    Examples
    --------
        >>> # Boolean network with only essential inputs
        >>> bn = random_network(N=10, n=2, allow_degenerate_functions=False)

        >>> # Classic NK-Kauffman network allowing degenerate rules
        >>> bn = random_network(N=10, n=3, allow_degenerate_functions=True)

        >>> # Fixed wiring: reuse an existing diagram but resample rules
        >>> bn0 = random_network(N=6, n=2)
        >>> bn  = random_network(I=bn)

        >>> # Exact canalizing depth k for all nodes
        >>> bn = random_network(N=8, n=3, depth=1, exact_depth=True)

        >>> # Nested canalizing update rules with specific layer structure (broadcast)
        >>> bn = random_network(N=5, n=3, layer_structure=[1,2])  # same for all nodes

        >>> # Parity rules
        >>> bn = random_network(N=7, n=2, parity=True)

        >>> # Poisson in-degrees (truncated), no self-regulation, request strong connectivity
        >>> bn = random_network(N=12, n=1.6, indegree_distribution='poisson',
        ...                     allow_self_loops=False, strongly_connected=True)

        >>> # Exact Hamming weights (broadcast)
        >>> bn = random_network(N=6, n=3, hamming_weight=4)

        >>> # To ensure strong connectivity, set allow_degenerate_functions=False
        >>> # and strongly_connected=True
        >>> bn = random_network(N,n,allow_degenerate_functions=False,strongly_connected=True)
    """
    rng = utils._coerce_rng(rng)
    if I is None and N is not None and n is not None:  # generate wiring diagram
        I = random_wiring_diagram(
            N,
            n,
            allow_self_loops=allow_self_loops,
            strongly_connected=strongly_connected,
            indegree_distribution=indegree_distribution,
            min_out_degree_one=min_out_degree_one,
            max_strong_connectivity_attempts=max_strong_connectivity_attempts,
            rng=rng,
        )

    elif I is not None:  # load wiring diagram
        assert isinstance(I, (list, np.ndarray, WiringDiagram, nx.DiGraph)), (
            "I must be an instance of WiringDiagram or a list or np.array of lists or np.arrays. Each inner list describes the regulators of node i (indexed by 0,1,...,len(I)-1)"
        )
        N = len(I)
        if isinstance(I, (list, np.ndarray)):
            for regulators in I:
                assert (
                    utils.is_list_or_array_of_ints(regulators)
                    and min(regulators) >= 0
                    and max(regulators) <= N - 1
                ), (
                    "Each element in I describes the regulators of a node (indexed by 0,1,...,len(I)-1)"
                )
            I = WiringDiagram(I)
        elif isinstance(I, nx.DiGraph):
            I = WiringDiagram.from_DiGraph( I )        
    else:
        raise AssertionError(
            "At a minimum, the wiring diagram I must be provided or the network size N and degree parameter n."
        )

    # Process the inputs, turn single inputs into vectors of length N

    # since layer_structure takes precedence over depth, this block needs to run before the depth block to ensure depth is a vector and not reset to a single value
    if layer_structure is None:
        layer_structure = [None] * N
    elif utils.is_list_or_array_of_ints(layer_structure):
        depth = sum(layer_structure)
        assert depth == 0 or (
            min(layer_structure) >= 1 and depth <= min(I.indegrees)
        ), (
            "The layer structure must be [] or a vector of positive integers with 0 <= depth = sum(layer_structure) <= N."
        )
        layer_structure = [layer_structure[:]] * N
    elif (
        np.all([utils.is_list_or_array_of_ints(el) for el in layer_structure])
        and len(layer_structure) == N
    ):
        for i, vector in enumerate(layer_structure):
            depth = sum(vector)
            assert depth == 0 or (min(vector) >= 1 and depth <= I.indegrees[i]), (
                "Ensure that layer_structure is an N-dimensional vector where each element represents a layer structure and is either [] or a vector of positive integers with 0 <= depth = sum(layer_structure[i]) <= n = indegrees[i]."
            )
    else:
        raise AssertionError(
            "Wrong input format for 'layer_structure'.\nIt must be a single vector (or N-dimensional vector of layer structures) where the sum of each element is between 0 and N."
        )

    if isinstance(depth, (int, np.integer)):
        assert depth >= 0, (
            "The canalizing depth must be an integer between 0 and min(indegrees) or an N-dimensional vector of integers must be provided to use different depths per function."
        )
        depth = [min(I.indegrees[i], depth) for i in range(N)]
    elif utils.is_list_or_array_of_ints(depth, required_length=N):
        depth = [min(I.indegrees[i], depth[i]) for i in range(N)]
        assert min(depth) >= 0, (
            "'depth' received a vector as input.\nTo use a user-defined vector, ensure that it is an N-dimensional vector where each element is a non-negative integer."
        )
    else:
        raise AssertionError(
            "Wrong input format for 'depth'.\nIt must be a single integer (or N-dimensional vector of integers) between 0 and N, specifying the minimal canalizing depth or exact canalizing depth (if exact_depth==True)."
        )

    if isinstance(bias, (float, np.floating)):
        bias = [bias] * N
    elif not utils.is_list_or_array_of_floats(bias, required_length=N):
        raise AssertionError(
            "Wrong input format for 'bias'.\nIt must be a single float (or N-dimensional vector of floats) in [0,1] , specifying the bias (probability of a 1) in the generation of the Boolean function."
        )

    if isinstance(absolute_bias, (float, np.floating)):
        absolute_bias = [absolute_bias] * N
    elif not utils.is_list_or_array_of_floats(absolute_bias, required_length=N):
        raise AssertionError(
            "Wrong input format for 'absolute_bias'.\nIt must be a single float (or N-dimensional vector of floats) in [0,1], specifying the absolute bias (divergence from the 'unbiased bias' of 0.5) in the generation of the Boolean function."
        )

    if hamming_weight == None:
        hamming_weight = [None] * N
    elif isinstance(hamming_weight, (int, np.integer)):
        hamming_weight = [hamming_weight] * N
    elif not utils.is_list_or_array_of_ints(hamming_weight, required_length=N):
        raise AssertionError(
            "Wrong input format for 'hamming_weight'.\nIf provided, it must be a single integer (or N-dimensional vector of integers) in {0,1,...,2^n}, specifying the number of 1s in the truth table of each Boolean function.\nIf exact_depth == True and depth==0, it must be in {2,3,...,2^n-2} because all functions with Hamming weight 0,1,2^n-1,2^n are canalizing."
        )

    # generate functions
    F = [
        random_function(
            n=I.indegrees[i],
            depth=depth[i],
            exact_depth=exact_depth,
            uniform_over_functions=uniform_over_functions,
            layer_structure=layer_structure[i],
            parity=parity,
            allow_degenerate_functions=allow_degenerate_functions,
            bias=bias[i],
            absolute_bias=absolute_bias[i],
            use_absolute_bias=use_absolute_bias,
            hamming_weight=hamming_weight[i],
            rng=rng,
        )
        for i in range(N)
    ]

    return BooleanNetwork(F, I)


def random_null_model(
    bn: BooleanNetwork,
    wiring_diagram: str = "fixed",
    preserve_bias: bool = True,
    preserve_canalizing_depth: bool = True,
    *,
    rng=None,
    **kwargs,
) -> BooleanNetwork:
    """
    Generate a randomized Boolean network (null model) from an existing
    Boolean network while preserving selected structural and dynamical
    properties.
        
    The returned network has the same number of nodes as ``bn``. Depending
    on the selected options, the wiring diagram and/or the Boolean update
    rules are randomized subject to specified invariants.
    
    Wiring diagram randomization
    ----------------------------
    The wiring diagram can be handled in one of three ways:
    
    - ``'fixed'`` (default):
      The original wiring diagram ``bn.I`` is reused unchanged.
    
    - ``'fixed_indegree'``:
      A new wiring diagram is sampled uniformly at random subject to
      preserving the in-degree of each node. This uses
      ``random_wiring_diagram`` with ``N = bn.N`` and ``n = bn.indegrees``.
    
    - ``'fixed_in_and_outdegree'``:
      The original wiring diagram is randomized via degree-preserving
      double-edge swaps using ``rewire_wiring_diagram``, preserving both
      in-degrees and out-degrees of all nodes.
    
    Rule randomization
    ------------------
    Independently of the wiring diagram, Boolean update rules are
    randomized for each node, optionally preserving properties of the
    original rules:
    
    - If ``preserve_bias`` is True, the exact Hamming weight (number of ones
      in the truth table) of each rule is preserved.
    
    - If ``preserve_canalizing_depth`` is True, the canalizing depth of each
      rule is preserved exactly.
    
    If both flags are True, both properties are preserved simultaneously.
    If neither flag is True, rules are regenerated subject only to
    non-degeneracy and the node’s in-degree.
    
    Parameters
    ----------
    bn : BooleanNetwork
        Source Boolean network.
    wiring_diagram : {'fixed', 'fixed_indegree', 'fixed_in_and_outdegree'}, optional
        Strategy for handling the wiring diagram. Default is ``'fixed'``.
    preserve_bias : bool, optional
        If True, preserve the exact Hamming weight of each Boolean rule.
        Default is True.
    preserve_canalizing_depth : bool, optional
        If True, preserve the exact canalizing depth of each Boolean rule.
        Default is True.
    rng : int, numpy.random.Generator, numpy.random.RandomState, random.Random, or None, optional
        Random number generator or seed specification. Passed to
        ``utils._coerce_rng``.
    kwargs
        Additional keyword arguments forwarded to the wiring-diagram
        randomization routine:
    
        - For ``wiring_diagram == 'fixed_indegree'``, forwarded to
          ``random_wiring_diagram`` (e.g., ``allow_self_loops``,
          ``strongly_connected``).
    
        - For ``wiring_diagram == 'fixed_in_and_outdegree'``, forwarded to
          ``rewire_wiring_diagram`` (e.g., ``average_swaps_per_edge``,
          ``allow_new_self_loops``, ``allow_self_loop_rewiring``).

    
    Returns
    -------
    BooleanNetwork
        A randomized Boolean network satisfying the selected invariants.
    
    Raises
    ------
    AssertionError
        If invalid options are provided.
    RuntimeError
        If wiring-diagram randomization fails (e.g., strong connectivity
        cannot be achieved within the allowed number of attempts).
    
    Notes
    -----
    This function generates null models by selectively preserving structural
    and dynamical properties of an existing Boolean network. It is intended
    for hypothesis testing and comparative studies rather than for uniform
    sampling over all networks satisfying the given constraints.
    
    Examples
    --------
    >>> # Most restrictive use case: Preserve both wiring and rule properties (default) 
    >>> bn_null = random_null_model(bn)
    
    >>> # Preserve in-degrees only and preserve rule bias
    >>> bn_null = random_null_model(
    ...     bn,
    ...     wiring_diagram='fixed_indegree',
    ...     preserve_bias=True,
    ...     preserve_canalizing_depth=False
    ... )
    
    >>> # Preserve both in- and out-degrees via rewiring
    >>> bn_null = random_null_model(
    ...     bn,
    ...     wiring_diagram='fixed_in_and_outdegree',
    ...     average_swaps_per_edge=15
    ... )
    """
    rng = utils._coerce_rng(rng)
    if wiring_diagram == "fixed":
        I = bn.I
    elif wiring_diagram == "fixed_indegree":
        I = random_wiring_diagram(N=bn.N, n=bn.indegrees, rng=rng, **kwargs)
    elif wiring_diagram == "fixed_in_and_outdegree":
        I = rewire_wiring_diagram(I=bn.I, **kwargs)
    else:
        raise AssertionError(
            "There are three choices for the wiring diagram: 1. 'fixed' (i.e., as in the provided BooleanNetwork), 2. 'fixed_indegree' (i.e., edges are shuffled but the indegree is preserved), 3. 'fixed_in_and_outdegree' (i.e., edges are shuffled but both the indegree and outdegree are preserved)."
        )

    dict_identity_nodes = bn.get_identity_nodes(as_dict=True)

    F = []
    for i, f in enumerate(bn.F):
        if dict_identity_nodes[i]:  # identity nodes don't change
            F.append(np.array([0, 1], dtype=int))
            continue
        if preserve_canalizing_depth:
            depth = f.get_canalizing_depth()
        if preserve_bias and preserve_canalizing_depth:
            core_function = f.properties["CoreFunction"]
            can_outputs = f.properties["CanalizedOutputs"]

            can_inputs = rng.choice(2, depth, replace=True)
            can_order = rng.choice(f.n, depth, replace=False)
            if f.n - depth == 0:
                core_function = np.array([1 - can_outputs[-1]], dtype=int)
            elif f.n - depth == 2:
                core_function = rng.choice(
                    [
                        np.array([0, 1, 1, 0], dtype=int),
                        np.array([1, 0, 0, 1], dtype=int),
                    ]
                )
            else:  # if f.n-depth>=3
                hamming_weight = sum(core_function)
                while True:
                    core_function = random_function_with_exact_hamming_weight(
                        f.n - depth, hamming_weight, rng=rng
                    )
                    if not core_function.is_canalizing():
                        if not core_function.is_degenerate():
                            break
            newf = np.full(2 ** bn.indegrees[i], -1, dtype=int)
            for j in range(depth):
                newf[
                    np.where(
                        np.bitwise_and(
                            newf == -1,
                            utils.get_left_side_of_truth_table(bn.indegrees[i])[
                                :, can_order[j]
                            ]
                            == can_inputs[j],
                        )
                    )[0]
                ] = can_outputs[j]
            newf[np.where(newf == -1)[0]] = core_function
            newf = BooleanFunction(newf)
        elif preserve_bias:  # and preserve_canalizing_depth==False
            hamming_weight = f.get_hamming_weight()
            newf = random_function_with_exact_hamming_weight(
                bn.indegrees[i], hamming_weight, rng=rng
            )
        elif preserve_canalizing_depth:
            newf = random_k_canalizing_function(
                n=bn.indegrees[i], k=depth, exact_depth=True, rng=rng
            )
        else:
            newf = random_non_degenerate_function(n=bn.indegrees[i], rng=rng)
        F.append(newf)
    return BooleanNetwork(F, I)
