#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections.abc import Sequence
import numpy as np

import utils

class MultistateFunction(object):
    def __init__(
            self,
            f : Sequence[int] | str,
            r : int,
            in_degree : Sequence[int],
            variables : Sequence[str] | None = None,
            name : str = ""):
        self.name = name
        if not isinstance(r, int):
            raise ValueError("Radix r must be an integer")
        if not r > 1:
            raise ValueError(f"Radix r has a minimum value of 2. Received {r}")
        if isinstance(f, str):
            f, r, self.variables = utils.from_from_expression(f)
        else:
            if not isinstance(f, (list, np.ndarray)):
                raise ValueError("f must be a list, numpy array, or interpretable string")
            if not len(f) > 0:
                raise ValueError("f cannot be empty")
            #_n = ???
            #if not len(f) == expected length:
            #    raise ValueError("f must be of size ???")
            if variables is None:
                self.variables = np.array([f"x{i}" for i in range(len(in_degree))])
            else:
                self.variables = np.asarray(variables, str)
                if self.variables.ndim != 1:
                    raise ValueError("variables must be a 1D array of strings")
                #if len(self.variables) != _n:
                #    raise ValueError(f"variables must have length {_n}, got {len(self.variables)}")
        self.n = len(self.variables)
        self.r = r
        self.in_rs = in_degree
        self.f = np.array(f, np.uint8)
        
        if not np.all(self.f < self.r):
            raise ValueError(f"f must contain only values 0 <= v < {self.r}")
#             _n = int(np.log2(len(f)))
#             if not abs(np.log2(len(f)) - _n) < 1e-9:
#                 raise ValueError("f must be of size 2^n, n >= 0")
            
#             if variables is None:
#                 self.variables = np.array([f"x{i}" for i in range(_n)])
#             else:
#                 self.variables = np.asarray(variables, dtype=str)
#                 if self.variables.ndim != 1:
#                     raise ValueError("variables must be a 1D array of strings")
#                 if len(self.variables) != _n:
#                     raise ValueError(
#                         f"variables must have length {_n}, got {len(self.variables)}"
#                     )
        
#         self.n = len(self.variables)
            
#         self.f = np.array(f, dtype=np.uint8)

#         if not np.all((self.f == 0) | (self.f == 1)):
#             raise ValueError("f must contain only the values 0 and 1.")
        
#         self.properties = {}


    ## Magic methods
    
    def __str__(self):
        """
        Return a human-readable string representation of the Boolean function.
        
        This method returns the underlying truth table as a NumPy array.
        """
        return f"{self.f}"
    
    def __repr__(self):
        """
        Return an unambiguous string representation of a MultistateFunction.
        
        For small functions (less than 6 variables), the full truth table is shown.
        For larger functions, only the number of inputs is displayed to
        avoid excessive output.
        """
        if self.n < 6:
            return f"{type(self).__name__}(f={self.f.tolist()})"
        return f"{type(self).__name__}(n={self.n})"
    
    def __len__(self):
        return len(self.f)
    
    def __getitem__(self, index):
        try:
            return int(self.f[index])
        except TypeError:
            return self.f[index]
    
    def __setitem__(self, index, value):
        self.f[index] = value

    # def __mul__(self, value)
    
    # def __rmul__(self, value):
    #     return self.__mul__(value)
    
    # def __and__(self, value)
    # def __or__(self, value)
    # def __xor__(self, value)
    # def __invert__(self)
    
#     def __call__(self, values: list[int] | tuple[int, ...] | np.ndarray):
#         """
#         Evaluate the Boolean function on a given input vector.
    
#         This method makes BooleanFunction instances callable and returns the
#         output value for a specified binary input configuration.
    
#         Parameters
#         ----------
#         values : list[int] | tuple[int, ...] | np.ndarray
#             Sequence of binary values (0 or 1) of length ``n``, where ``n`` is
#             the number of input variables of the Boolean function.
    
#         Returns
#         -------
#         int
#             Output value of the Boolean function (0 or 1) for the specified input.
    
#         Raises
#         ------
#         ValueError
#             If the input length does not match ``n`` or if non-binary values are
#             provided.
    
#         Examples
#         --------
#         >>> f = BooleanFunction("x1 | (x2 & x3)")
#         >>> f([1, 0, 1])
#         1
#         >>> f([0, 1, 0])
#         0
#         """
#         if not len(values) == self.n:
#             raise ValueError(f"The argument must be of length {self.n}.")
#         if not set(values) <= {0, 1}:
#             raise ValueError("Binary values required.")
#         return self.f[utils.bin2dec(values)].item()
#     ## Conversions:
    

#     def to_polynomial(self) -> str:
#         """
#         Convert the Boolean function to a polynomial representation.
    
#         This method returns a polynomial representation of the Boolean function
#         in non-reduced disjunctive normal form (DNF).
    
#         Returns
#         -------
#         str
#             Polynomial representation of the Boolean function in non-reduced DNF.
#         """
#         return utils.bool_to_poly(self.f, variables=self.variables)


#     def to_truth_table(
#         self,
#         return_output: bool = True,
#         filename: str | None = None,
#     ) -> pd.DataFrame | None:
#         """
#         Return or save the full truth table of the Boolean function.
    
#         The truth table is represented as a pandas DataFrame in which each row
#         corresponds to an input configuration and the final column contains the
#         output value of the Boolean function.
    
#         Parameters
#         ----------
#         return_output : bool, optional
#             Whether to return the truth table as a DataFrame. If ``False``, the
#             truth table is only written to file when ``filename`` is provided.
#             Default is ``True``.
#         filename : str or None, optional
#             File name (including extension) to which the truth table is saved.
#             Supported formats are ``'csv'``, ``'xls'``, and ``'xlsx'``. If
#             provided, the truth table is automatically written to disk.
    
#         Returns
#         -------
#         pandas.DataFrame or None
#             The full truth table if ``return_output=True``; otherwise ``None``.
    
#         Notes
#         -----
#         The column names correspond to the input variable names followed by the
#         function name if provided, or ``'f'`` otherwise. When saving to a file,
#         the output format is determined by the file extension.
    
#         Examples
#         --------
#         >>> f = BooleanFunction("(x1 & ~x2) | x3")
#         >>> f.to_truth_table()
#            x1  x2  x3  f
#         0   0   0   0  0
#         1   0   0   1  1
#         2   0   1   0  0
#         3   0   1   1  1
#         4   1   0   0  1
#         5   1   0   1  1
#         6   1   1   0  0
#         7   1   1   1  1
#         """
#         columns = np.append(self.variables, self.name if self.name != "" else "f")
#         truth_table = pd.DataFrame(
#             np.c_[utils.get_left_side_of_truth_table(self.n), self.f],
#             columns=columns,
#         )
    
#         if filename is not None:
#             ending = filename.split(".")[-1]
#             if ending not in {"csv", "xls", "xlsx"}:
#                 raise ValueError("filename must end in 'csv', 'xls', or 'xlsx'")
#             if ending == "csv":
#                 truth_table.to_csv(filename)
#             else:
#                 truth_table.to_excel(filename)
    
#         if return_output:
#             return truth_table
#         return None


#     def to_logical(
#         self,
#         and_op: str = "&",
#         or_op: str = "|",
#         not_op: str = "!",
#         minimize_expression: bool = True,
#     ) -> str:
#         """
#         Convert the Boolean function to a logical expression.
    
#         This method converts the Boolean function from its truth table
#         representation to a logical expression. If the PyEDA package is
#         available, the expression can optionally be minimized using the
#         Espresso algorithm. Otherwise, a non-minimized expression is
#         generated as a fallback.
    
#         Parameters
#         ----------
#         and_op : str, optional
#             String used to represent the logical AND operator. Default is ``"&"``.
#         or_op : str, optional
#             String used to represent the logical OR operator. Default is ``"|"``.
#         not_op : str, optional
#             String used to represent the logical NOT operator. Default is ``"!"``.
#         minimize_expression : bool, optional
#             Whether to minimize the logical expression using the Espresso
#             algorithm (via PyEDA). If ``False``, the expression is returned
#             in non-minimized disjunctive normal form. Default is ``True``.
    
#         Returns
#         -------
#         str
#             Logical expression representing the Boolean function.
    
#         Notes
#         -----
#         If the PyEDA package is not installed, the method falls back to a
#         non-minimized expression derived from the polynomial representation.
#         """
#         try:
#             from pyeda.inter import exprvar, Or, And, espresso_exprs
#             from pyeda.boolalg.expr import OrOp, AndOp, NotOp, Complement
#             __LOADED_PyEDA__ = True
#         except ModuleNotFoundError:
#             __LOADED_PyEDA__ = False
            
#         if __LOADED_PyEDA__:
#             variables = [exprvar(str(var)) for var in self.variables]
#             minterms = [i for i in range(2**self.n) if self.f[i]]
    
#             terms = []
#             for m in minterms:
#                 bits = [
#                     variables[i]
#                     if (m >> (self.n - 1 - i)) & 1
#                     else ~variables[i]
#                     for i in range(self.n)
#                 ]
#                 terms.append(And(*bits))
    
#             func_expr = Or(*terms).to_dnf()
    
#             if func_expr.is_zero():
#                 return "0"
    
#             if minimize_expression:
#                 func_expr, = espresso_exprs(func_expr)
    
#             def __pyeda_to_string__(e):
#                 if isinstance(e, OrOp):
#                     return "(" + f"){or_op}(".join(
#                         __pyeda_to_string__(arg) for arg in e.xs
#                     ) + ")"
#                 elif isinstance(e, AndOp):
#                     return and_op.join(__pyeda_to_string__(arg) for arg in e.xs)
#                 elif isinstance(e, NotOp):
#                     return f"{not_op}({__pyeda_to_string__(e.x)})"
#                 elif isinstance(e, Complement):
#                     return f"({not_op}{str(e)[1:]})"
#                 return str(e)
    
#             return __pyeda_to_string__(func_expr)
#         else:
#             # Fallback without PyEDA
#             return (
#                 self.to_polynomial()
#                 .replace(" * ", and_op)
#                 .replace(" + ", or_op)
#                 .replace("1 - ", not_op)
#             )

    
#     def summary(self, *, as_dict: bool = False, compute_all: bool = False):
#         """
#         Return a concise summary of the Boolean function.
    
#         The summary includes basic structural and statistical properties of the
#         Boolean function and, optionally, additional properties that may require
#         nontrivial computation.
    
#         Parameters
#         ----------
#         as_dict : bool, optional
#             If ``True``, return the summary as a dictionary. If ``False`` (default),
#             return a formatted string.
    
#         compute_all : bool, optional
#             If ``True``, additional properties are computed and included in the
#             summary. These computations may be expensive. If ``False`` (default),
#             only already available properties are included.
    
#         Returns
#         -------
#         str or dict
#             Summary of the Boolean function, either as a formatted string or as
#             a dictionary depending on the value of ``as_dict``.
#         """
    
#         ones = sum(self.f)
#         bias = ones / 2**self.n
#         abs_bias = 2 * abs(bias-0.5)

#         summary = {
#             "Number of variables": self.n,
#             "Hamming Weight": int(ones),
#             "Bias": float(bias),
#             "Absolute bias": float(abs_bias),
#             "Variables": self.variables
#         }
    
#         # Optional properties (only if already computed / cached)
#         if compute_all:
#             self.get_layer_structure()
#             self.get_type_of_inputs()
            
#         summary.update(self.properties)
                
#         if as_dict:
#             return summary
    
#         # Pretty formatting
#         lines = [
#             "BooleanFunction summary",
#             "-" * 40,
#             f"Number of variables:        {summary['Number of variables']}",
#             f"Hamming Weight:             {summary['Hamming Weight']}",
#             f"Bias:                       {summary['Bias']:.3f}",
#             f"Absolute bias:              {summary['Absolute bias']:.3f}",
#             f"Variables:                  {summary['Variables']}"
#         ]
    
#         for key in summary:
#             if key not in {
#                 "Number of variables",
#                 "Hamming Weight",
#                 "Bias",
#                 "Absolute bias",
#                 "Variables"
#             }:
#                 lines.append(f"{key}:"+(" "*(27-len(key)))+f"{summary[key]}")

#         return "\n".join(lines)

    
#     def is_constant(self) -> bool:
#         """
#         Check whether the Boolean function is constant.
    
#         A Boolean function is constant if all entries of its truth table are
#         identical (all 0 or all 1).
    
#         Returns
#         -------
#         bool
#             ``True`` if the Boolean function is constant, ``False`` otherwise.
#         """
#         return bool(np.all(self.f == self.f[0]))
        
#     def is_degenerate(self, use_numba: bool = True) -> bool:
#         """
#         Determine whether the Boolean function is degenerate.
    
#         A Boolean function is degenerate if it contains at least one
#         non-essential variable, i.e., a variable on which the function's
#         output does not depend.
    
#         Parameters
#         ----------
#         use_numba : bool, optional
#             Whether to use Numba-accelerated computation when available.
#             Default is ``True``.
    
#         Returns
#         -------
#         bool
#             ``True`` if the Boolean function contains at least one
#             non-essential variable, ``False`` if all variables are essential.
#         """
#         if __LOADED_NUMBA__ and use_numba:
#             f = np.asarray(self.f, dtype=np.uint8)
#             return bool(_is_degenerate_numba(f, self.n))
#         else:
#             for i in range(self.n):
#                 dummy_add = 2 ** (self.n - 1 - i)
#                 dummy = np.arange(2**self.n) % (2 ** (self.n - i)) // dummy_add
#                 depends_on_i = False
#                 for j in range(2**self.n):
#                     if dummy[j] == 1:
#                         continue
#                     else:
#                         if self.f[j] != self.f[j + dummy_add]:
#                             depends_on_i = True
#                             break
#                 if not depends_on_i:
#                     return True
#             return False
        
        
#     def is_monotonic(self) -> bool:
#         """
#         Determine whether the Boolean function is monotonic.
    
#         A Boolean function is monotonic if it is monotonic in each variable,
#         i.e., for every variable the function is either non-decreasing or
#         non-increasing with respect to that variable.
    
#         Returns
#         -------
#         bool
#             ``True`` if the Boolean function is monotonic, ``False`` otherwise.
#         """
#         return "conditional" not in self.get_type_of_inputs()
    
    
#     def is_canalizing(self) -> bool:
#         """
#         Determine whether the Boolean function is canalizing.
    
#         A Boolean function is canalizing if there exists at least one variable
#         and a value in ``{0,1}`` such that fixing that variable to the given
#         value forces the output of the function to be constant.
    
#         Returns
#         -------
#         bool
#             ``True`` if the Boolean function is canalizing, ``False`` otherwise.
#         """
#         indices = np.arange(2**self.n, dtype=np.uint32)
    
#         for i in range(self.n):
#             mask = 1 << i
#             bit_is_0 = (indices & mask) == 0
#             bit_is_1 = ~bit_is_0
    
#             f0 = self.f[bit_is_0]
#             f1 = self.f[bit_is_1]
    
#             if np.all(f0 == f0[0]) or np.all(f1 == f1[0]):
#                 return True
    
#         return False

#     def is_k_canalizing(self, k: int) -> bool:
#         """
#         Determine whether the Boolean function is k-canalizing.
    
#         A Boolean function is k-canalizing if it has a sequence of at least
#         ``k`` canalizing variables. After fixing a canalizing variable to its
#         canalizing value, the resulting subfunction must itself be
#         (k−1)-canalizing, recursively.
    
#         Parameters
#         ----------
#         k : int
#             Desired canalizing depth, with ``0 <= k <= n``. Every Boolean
#             function is trivially 0-canalizing.
    
#         Returns
#         -------
#         bool
#             ``True`` if the Boolean function is k-canalizing, ``False`` otherwise.
    
#         Notes
#         -----
#         This method has exponential time complexity in ``n`` and is intended for
#         small Boolean functions.
    
#         References
#         ----------
#         He, Q., & Macauley, M. (2016).
#             Stratification and enumeration of Boolean functions by canalizing depth.
#             *Physica D: Nonlinear Phenomena*, 314, 1–8.
    
#         Dimitrova, E., Stigler, B., Kadelka, C., & Murrugarra, D. (2022).
#             Revealing the canalizing structure of Boolean functions:
#             Algorithms and applications.
#             *Automatica*, 146, 110630.
#         """
#         if k > self.n:
#             return False
#         if k == 0:
#             return True
#         if np.all(self.f == self.f[0]):
#             return False
    
#         indices = np.arange(2**self.n, dtype=np.uint32)
    
#         for i in range(self.n):
#             mask = 1 << i
#             bit_is_0 = (indices & mask) == 0
#             bit_is_1 = ~bit_is_0
    
#             f0, f1 = self.f[bit_is_0], self.f[bit_is_1]
    
#             if np.all(f0 == f0[0]):
#                 return True if k == 1 else BooleanFunction(f1).is_k_canalizing(k - 1)
    
#             if np.all(f1 == f1[0]):
#                 return True if k == 1 else BooleanFunction(f0).is_k_canalizing(k - 1)
    
#         return False

#     def is_kset_canalizing(self, k: int) -> bool:
#         """
#         Determine whether the Boolean function is k-set canalizing.
    
#         A Boolean function is k-set canalizing if there exists a set of ``k``
#         variables such that fixing these variables to specific values forces
#         the output of the function, regardless of the remaining ``n - k``
#         variables.
    
#         Parameters
#         ----------
#         k : int
#             Size of the variable set, with ``0 <= k <= n``.
    
#         Returns
#         -------
#         bool
#             ``True`` if the Boolean function is k-set canalizing, ``False`` otherwise.
    
#         Notes
#         -----
#         This method has exponential time complexity in ``n`` and is intended for
#         small Boolean functions.
    
#         References
#         ----------
#         Kadelka, C., Keilty, B., & Laubenbacher, R. (2023).
#             Collectively canalizing Boolean functions.
#             *Advances in Applied Mathematics*, 145, 102475.
#         """
#         return self.get_kset_canalizing_proportion(k) > 0


#     ## Methods with non-binary output

#     def get_hamming_weight(self) -> int:
#         """
#         Compute the Hamming weight of the Boolean function.
    
#         The Hamming weight is the number of input states for which the function
#         evaluates to ``1`` (i.e., the number of ones in the truth table).
    
#         Returns
#         -------
#         int
#             The Hamming weight of the Boolean function.
#         """
#         return int(self.f.sum())
        
    
#     def get_essential_variables(self) -> list:
#         """
#         Determine the essential variables of the Boolean function.
    
#         A variable ``x_i`` is essential if there exists at least one assignment of
#         the remaining variables such that flipping ``x_i`` changes the output of
#         the function.
    
#         Returns
#         -------
#         list[int]
#             Indices of all essential variables. If the truth table is empty, returns
#             an empty list.
#         """
#         if len(self.f) == 0:
#             return []
#         essential_variables = list(range(self.n))
#         for i in range(self.n):
#             dummy_add = (2**(self.n-1-i))
#             dummy = np.arange(2**self.n) % (2**(self.n-i)) // dummy_add
#             depends_on_i = False
#             for j in range(2**self.n):
#                 if dummy[j] == 1:
#                     continue
#                 else:
#                     if self.f[j] != self.f[j + dummy_add]:
#                         depends_on_i = True
#                         break
#             if depends_on_i == False:
#                 essential_variables.remove(i)
#         return essential_variables 

#     def get_number_of_essential_variables(self) -> int:
#         """
#         Count the number of essential variables of the Boolean function.
    
#         Returns
#         -------
#         int
#             The number of essential variables.
#         """
#         return len(self.get_essential_variables())
    
    
#     def get_type_of_inputs(self) -> np.ndarray:
#         """
#         Classify each input variable of the Boolean function.
    
#         Each variable is classified as one of:
    
#         - ``'non-essential'``: flipping the variable never changes the output
#         - ``'positive'``: flipping the variable from 0 to 1 never decreases the output
#         - ``'negative'``: flipping the variable from 0 to 1 never increases the output
#         - ``'conditional'``: flipping the variable can both increase and decrease the output
    
#         The result is cached in ``self.properties['InputTypes']``.
    
#         Returns
#         -------
#         np.ndarray
#             Array of shape ``(n,)`` with dtype ``str`` giving the type of each input
#             variable.
#         """

#         if 'InputTypes' in self.properties:
#             return self.properties['InputTypes']
    
#         f = np.asarray(self.f, dtype=np.int8)
#         n = self.n
    
#         types = np.empty(n, dtype=object)
    
#         # Compute all pairwise differences for each bit position simultaneously # Each variable toggles every 2**i entries in the truth table. 
#         for i in range(n): 
#             period = 2 ** (i + 1) 
#             half = period // 2 
            
#             # Vectorized reshape pattern: consecutive blocks of 0...1 transitions 
#             f_reshaped = f.reshape(-1, period) 
#             diff = f_reshaped[:, half:] - f_reshaped[:, :half] 
#             min_diff = diff.min() 
#             max_diff = diff.max() 
#             if min_diff == 0 and max_diff == 0: 
#                 types[i] = 'non-essential' 
#             elif min_diff == -1 and max_diff == 1: 
#                 types[i] = 'conditional' 
#             elif min_diff >= 0 and max_diff == 1: 
#                 types[i] = 'positive' 
#             elif min_diff == -1 and max_diff <= 0: 
#                 types[i] = 'negative'
                
#         types = np.array(types, dtype=str)
#         self.properties['InputTypes'] = types
#         return types[::-1] #flip because of BoolForge logic ordering
        
    
#     def get_symmetry_groups(self) -> list[list[int]]:
#         """
#         Identify symmetry groups of input variables.
    
#         Two variables belong to the same symmetry group if swapping their values
#         leaves the Boolean function invariant for all assignments of the remaining
#         variables.
    
#         Returns
#         -------
#         list[list[int]]
#             A list of symmetry groups, where each group is given by a list of
#             variable indices.
            
#         Notes
#         -----
#         This method has exponential time complexity in ``n`` and is intended for
#         small Boolean functions.
#         """
#         symmetry_groups = []
#         left_to_check = np.ones(self.n)
#         for i in range(self.n):
#             if left_to_check[i] == 0:
#                 continue
#             else:
#                 symmetry_groups.append([i])
#                 left_to_check[i] = 0
#             for j in range(i + 1, self.n):
#                 diff = sum(2**np.arange(self.n - i - 2, self.n - j - 2, -1))
#                 for ii, x in enumerate(utils.get_left_side_of_truth_table(self.n)):
#                     if x[i] != x[j] and x[i] == 0 and self.f[ii] != self.f[ii + diff]:
#                         break
#                 else:
#                     left_to_check[j] = 0
#                     symmetry_groups[-1].append(j)
#         return symmetry_groups
    
#     def get_absolute_bias(self) -> float:
#         """
#         Compute the absolute bias of the Boolean function.
    
#         The absolute bias is defined as
    
#         ``| (H / 2**(n-1)) - 1 |``,
    
#         where ``H`` is the Hamming weight of the function. It measures how far the
#         output distribution deviates from being perfectly balanced.
    
#         Returns
#         -------
#         float
#             The absolute bias of the Boolean function.
#         """
#         return float(abs(self.get_hamming_weight() * 1.0 / 2**(self.n - 1) - 1))


#     def get_activities(
#         self,
#         nsim: int = 10000,
#         exact: bool = False,
#         *,
#         rng=None
#     ) -> np.ndarray:
#         """
#         Compute the activities of all input variables.
    
#         The activity of a variable is the probability that flipping this variable
#         (while keeping all others fixed) changes the output of the Boolean function.
    
#         Activities can be computed exactly by enumerating all ``2**n`` input states
#         or estimated via Monte Carlo sampling.
    
#         Parameters
#         ----------
#         nsim : int, optional
#             Number of random samples used when ``exact=False`` (default: 10000).
#         exact : bool, optional
#             If ``True``, compute activities exactly by enumerating all input states.
#             If ``False``, estimate activities via sampling (default: ``False``).
#         rng : None or numpy.random.Generator, optional
#             Random number generator passed to ``utils._coerce_rng``.
    
#         Returns
#         -------
#         np.ndarray
#             Array of shape ``(n,)`` containing the activities of all variables.
#         """     
#         size_state_space = 2**self.n
#         activities = np.zeros(self.n,dtype=np.float64)
#         if exact:
#             # Compute all integer representations of inputs (0 .. 2^n - 1)
#             X = np.arange(size_state_space, dtype=np.uint32)
        
#             # For each bit position i, flipping that bit corresponds to XOR with (1 << self.n-1-i)
#             for i in range(self.n):
#                 flipped = X ^ (1 << self.n-1-i) 
#                 activities[i] = np.count_nonzero(self.f != self.f[flipped])
#             return activities / size_state_space
#         else:
#             rng = utils._coerce_rng(rng)

#             random_states = rng.integers(0,size_state_space,nsim) #
#             for i in range(self.n):
#                 flipped_random_states = random_states ^ (1 << self.n-1-i) 
#                 activities[i] = np.count_nonzero(self.f[random_states] != self.f[flipped_random_states])
#             return activities / nsim
    
    
#     def get_average_sensitivity(
#         self,
#         nsim: int = 10000,
#         exact: bool = False,
#         normalized: bool = True,
#         *,
#         rng=None
#     ) -> float:
#         """
#         Compute the average sensitivity of the Boolean function.
    
#         The (unnormalized) average sensitivity equals the sum of the activities of
#         all variables. If ``normalized=True``, the result is divided by ``n``.
    
#         The sensitivity can be computed exactly by enumerating all input states or
#         estimated via Monte Carlo sampling.
    
#         Parameters
#         ----------
#         nsim : int, optional
#             Number of random samples used when ``exact=False`` (default: 10000).
#         exact : bool, optional
#             If ``True``, compute the exact activities by enumerating all input states.
#             If ``False``, estimate them via sampling (default: ``False``).
#         normalized : bool, optional
#             If ``True`` (default), return the average sensitivity divided by ``n``.
#             If ``False``, return the sum of activities.
#         rng : None or numpy.random.Generator, optional
#             Random number generator passed to ``utils._coerce_rng``.
    
#         Returns
#         -------
#         float
#             The (optionally normalized) average sensitivity of the Boolean function.
#         """      
#         activities = self.get_activities(nsim,exact,rng=rng)
#         s = sum(activities)
#         if normalized:
#             return float(s / self.n)
#         else:
#             return float(s)
    

#     def _get_layer_structure(
#         self,
#         can_inputs,
#         can_outputs,
#         can_order,
#         variables,
#         depth,
#         number_layers
#     ):
#         """
#         Internal recursive routine for computing the canalizing layer structure.
    
#         This method identifies all canalizing variables at the current recursion
#         level using bitwise operations, removes them simultaneously, and recurses
#         on the resulting core function.
    
#         Parameters
#         ----------
#         can_inputs : np.ndarray
#             Accumulated canalizing input values.
#         can_outputs : np.ndarray
#             Accumulated canalized output values.
#         can_order : np.ndarray
#             Accumulated order of canalizing variables.
#         variables : list[int]
#             Indices of variables remaining in the current subfunction.
#         depth : int
#             Current canalizing depth.
#         number_layers : int
#             Current number of identified canalizing layers.
    
#         Returns
#         -------
#         tuple
#             A tuple containing the updated canalizing depth, number of layers,
#             canalizing inputs, canalized outputs, core Boolean function, and
#             canalizing variable order.
#         """

#         # base cases
#         if np.all(self.f == self.f[0]):
#             # recursion ends when function becomes constant
#             return depth, number_layers, can_inputs, can_outputs, self, can_order
    
#         if not variables:
#             variables = list(range(self.n))
#         elif isinstance(variables, np.ndarray):
#             variables = variables.tolist()
    
#         indices = np.arange(2**self.n, dtype=np.uint32)
    
#         # candidate canalizing variables (x_i, a)
#         new_canalizing_vars = []
#         new_can_inputs = []
#         new_can_outputs = []
#         new_f = None
    
#         for i in range(self.n):
#             mask = 1 << (self.n-1-i)
#             bit0 = (indices & mask) == 0
#             bit1 = ~bit0
#             f0, f1 = self.f[bit0], self.f[bit1]
    
#             # check both possible canalizing directions
#             if np.all(f0 == f0[0]):
#                 new_canalizing_vars.append(variables[i])
#                 new_can_inputs.append(0)
#                 new_can_outputs.append(int(f0[0]))
#             elif np.all(f1 == f1[0]):
#                 new_canalizing_vars.append(variables[i])
#                 new_can_inputs.append(1)
#                 new_can_outputs.append(int(f1[0]))
    
#         if not new_canalizing_vars:
#             # non-canalizing core function
#             return depth, number_layers, can_inputs, can_outputs, self, can_order
    
#         # reduce variable list (remove canalizing ones)
#         indices_new_canalizing_vars = [i for i,v in enumerate(variables) if v in new_canalizing_vars]
#         remaining_vars = [v for v in variables if v not in new_canalizing_vars]
    
#         # build the restricted subfunction (“core function”)
#         # start with all indices, then keep those where none of the canalizing
#         # variables take their canalizing inputs
#         mask_keep = np.ones(2**self.n, dtype=bool)
#         for var, val in zip(indices_new_canalizing_vars, new_can_inputs):
#             bitmask = (indices >> (self.n-1-var)) & 1
#             mask_keep &= (bitmask != val)
#         new_f = self.f[mask_keep]
    

#         # recurse on reduced function
#         new_bf = self.__class__(list(new_f))
#         return new_bf._get_layer_structure(
#             np.append(can_inputs, new_can_inputs),
#             np.append(can_outputs, new_can_outputs),
#             np.append(can_order, new_canalizing_vars),
#             remaining_vars,
#             depth + len(new_canalizing_vars),
#             number_layers + 1,
#         )


#     def get_layer_structure(self) -> dict:
#         """
#         Determine the canalizing layer structure of a Boolean function.
        
#         This method decomposes a Boolean function into its canalizing layers
#         (standard monomial form) by recursively identifying and removing
#         canalizing variables. All variables that canalize the function at the
#         same recursion step form one canalizing layer and are removed
#         simultaneously.
        
#         The decomposition yields the canalizing depth, the number of canalizing
#         layers, the canalizing inputs and outputs, the order of canalizing
#         variables, and the remaining non-canalizing core function.
        
#         Returns
#         -------
#         dict
#             Dictionary containing the canalizing layer structure with the
#             following entries:
        
#             - ``CanalizingDepth`` : int  
#               Total number of canalizing variables.
        
#             - ``NumberOfLayers`` : int  
#               Number of distinct canalizing layers.
        
#             - ``CanalizingInputs`` : np.ndarray  
#               Canalizing input value for each canalizing variable.
        
#             - ``CanalizedOutputs`` : np.ndarray  
#               Output value forced by each canalizing variable.
        
#             - ``CoreFunction`` : BooleanFunction  
#               Core Boolean function obtained after removing all canalizing
#               variables.
        
#             - ``OrderOfCanalizingVariables`` : np.ndarray  
#               Order in which canalizing variables are identified.
        
#             - ``LayerStructure`` : np.ndarray  
#               Number of canalizing variables in each layer.
        
#         Notes
#         -----
#         The result is cached in ``self.properties`` and recomputed only if the
#         canalizing structure has not been computed previously.
        
#         Notes
#         -----
#         This method has exponential time complexity in ``n`` and is intended for
#         smaller Boolean functions.
        
#         References
#         ----------
#         He, Q., & Macauley, M. (2016).
#             Stratification and enumeration of Boolean functions by canalizing depth.
#             *Physica D: Nonlinear Phenomena*, 314, 1–8.
        
#         Dimitrova, E., Stigler, B., Kadelka, C., & Murrugarra, D. (2022).
#             Revealing the canalizing structure of Boolean functions:
#             Algorithms and applications.
#             *Automatica*, 146, 110630.
#         """
#         if "CanalizingDepth" not in self.properties:
#             dummy = dict(zip(["CanalizingDepth", "NumberOfLayers", "CanalizingInputs", "CanalizedOutputs", "CoreFunction", "OrderOfCanalizingVariables"],
#                                             self._get_layer_structure(can_inputs=np.array([], dtype=int), can_outputs=np.array([], dtype=int),
#                                                                       can_order=np.array([], dtype=int), variables=[], depth=0, number_layers=0)))
#             dummy.update({'LayerStructure': get_layer_structure_from_canalized_outputs(dummy["CanalizedOutputs"])})
#             self.properties.update(dummy)
#             return dummy
#         else:
#             return {key: self.properties[key] for key in ["CanalizingDepth", "NumberOfLayers", "CanalizingInputs", "CanalizedOutputs", "CoreFunction", "OrderOfCanalizingVariables",'LayerStructure']}


#     def get_canalizing_depth(self) -> int:
#         """
#         Return the canalizing depth of the Boolean function.
    
#         The canalizing depth is the total number of canalizing variables identified
#         in the canalizing layer decomposition.
    
#         Returns
#         -------
#         int
#             Canalizing depth of the Boolean function.
#         """
#         if "CanalizingDepth" not in self.properties:
#             self.get_layer_structure()
#         return self.properties["CanalizingDepth"]

    
#     def get_kset_canalizing_proportion(self, k : int) -> float:
#         """
#         Compute the proportion of k-set canalizing input sets.
    
#         For a given ``k``, this method computes the probability that a randomly
#         chosen set of ``k`` variables canalizes the function, i.e., fixing those
#         variables to some values forces the output regardless of the remaining
#         variables.
    
#         Parameters
#         ----------
#         k : int
#             Size of the variable set, with ``0 <= k <= n``.
    
#         Returns
#         -------
#         float
#             Proportion of k-set canalizing input sets.
            
#         Notes
#         -----
#         This method has exponential time complexity in ``n`` and is intended for
#         small Boolean functions.
    
#         References
#         ----------
#         Kadelka, C., Keilty, B., & Laubenbacher, R. (2023).
#             Collectively canalizing Boolean functions.
#             *Advances in Applied Mathematics*, 145, 102475.
#         """
#         if not (type(k)==int and 0<=k<=self.n):
#             raise ValueError("k must be an integer and satisfy 0 <= k <= degree n")
        
#         # trivial case
#         if k == 0:
#             return float(self.is_constant())
        
#         # precompute binary representation of all inputs
#         #indices = np.arange(2**self.n, dtype=np.uint32)
#         #bits = ((indices[:, None] >> np.arange(self.n)) & 1).astype(np.uint8)  # shape (2**n, n)
#         left_side_of_truth_table = utils.get_left_side_of_truth_table(self.n)
        
#         total_tests = 0
#         canalizing_hits = 0
    
#         # iterate over variable subsets of size k
#         for subset in itertools.combinations(range(self.n), k):
#             Xsub = left_side_of_truth_table[:, subset]  # shape (2**n, k)
#             # For each possible assignment to this subset
#             for assignment in itertools.product([0, 1], repeat=k):
#                 mask = np.all(Xsub == assignment, axis=1)
#                 if not np.any(mask):
#                     continue
#                 # If all outputs equal when these vars are fixed → canalizing
#                 f_sub = self.f[mask]
#                 if np.all(f_sub == f_sub[0]):
#                     canalizing_hits += 1
#                 total_tests += 1
    
#         return canalizing_hits / total_tests if total_tests > 0 else 0.0


#     def get_kset_canalizing_proportion_of_variables(self, k : int) -> float:
#         """
#         Compute the proportion of k-set canalizing input sets per variable.
    
#         For a given ``k``, this method computes, for each variable, the proportion
#         of k-variable input sets containing that variable which canalize the
#         Boolean function.
    
#         Parameters
#         ----------
#         k : int
#             Size of the variable set, with ``0 <= k <= n``.
    
#         Returns
#         -------
#         np.ndarray
#             Array of length ``n`` giving the proportion of k-set canalizing input
#             sets containing each variable.
            
#         Notes
#         -----
#         This method has exponential time complexity in ``n`` and is intended for
#         small Boolean functions.
    
#         References
#         ----------
#         Kadelka, C., Keilty, B., & Laubenbacher, R. (2023).
#             Collectively canalizing Boolean functions.
#             *Advances in Applied Mathematics*, 145, 102475.
#         """
#         if not (type(k)==int and 0<=k<=self.n):
#             raise ValueError("k must be an integer and satisfy 0 <= k <= degree n")
        
#         # trivial case
#         if k == 0:
#             return float(self.is_constant())
        
#         # precompute binary representation of all inputs
#         #indices = np.arange(2**self.n, dtype=np.uint32)
#         #bits = ((indices[:, None] >> np.arange(self.n)) & 1).astype(np.uint8)  # shape (2**n, n)
#         left_side_of_truth_table = utils.get_left_side_of_truth_table(self.n)
        
#         canalizing_hits = np.zeros(self.n,dtype=np.float64)
        
#         # iterate over variable subsets of size k
#         for subset in itertools.combinations(range(self.n), k):
#             Xsub = left_side_of_truth_table[:, subset]  # shape (2**n, k)
#             subset = np.array(subset)
#             # For each possible assignment to this subset
#             for assignment in itertools.product([0, 1], repeat=k):
#                 mask = np.all(Xsub == assignment, axis=1)
#                 if not np.any(mask):
#                     continue
#                 # If all outputs equal when these vars are fixed → canalizing
#                 f_sub = self.f[mask]
#                 if np.all(f_sub == f_sub[0]):
#                     canalizing_hits[subset] += 1
    
#         return canalizing_hits / (k/self.n * math.comb(self.n,k) * 2**k)


#     def get_canalizing_strength(self) -> tuple:
#         """
#         Compute the canalizing strength of the Boolean function.
    
#         The canalizing strength is defined as a weighted average of the proportions
#         of k-set canalizing inputs for ``k = 1, ..., n-1``. It equals 0 for minimally
#         canalizing functions (e.g., parity functions) and 1 for maximally canalizing
#         functions (e.g., nested canalizing functions with a single layer).
    
#         Returns
#         -------
#         float
#             Canalizing strength of the Boolean function.

#         Notes
#         -----
#         This method has exponential time complexity in ``n`` and is intended for
#         small Boolean functions.

#         References
#         ----------
#         Kadelka, C., Keilty, B., & Laubenbacher, R. (2023).
#             Collectively canalizing Boolean functions.
#             *Advances in Applied Mathematics*, 145, 102475.
#         """

#         if self.n==1:
#             warnings.warn(
#                 "Canalizing strength is only defined for Boolean functions with n > 1 inputs. "
#                 "Returning 1 for n == 1.",
#                 RuntimeWarning
#             )
#             return 1.0
#         res = []
#         for k in range(1, self.n):
#             res.append(self.get_kset_canalizing_proportion(k))
#         return float(np.mean(np.multiply(res, 2**np.arange(1, self.n) / (2**np.arange(1, self.n) - 1))))


#     def get_canalizing_strength_of_variables(self) -> np.ndarray:
#         """
#         Compute the canalizing strength of each variable.
    
#         The canalizing strength of a variable is defined as a weighted average of
#         the proportions of k-set canalizing inputs containing that variable for
#         ``k = 1, ..., n-1``.
    
#         Notes
#         -----
#         This method has exponential time complexity in ``n`` and is intended for
#         small Boolean functions.

#         Returns
#         -------
#         np.ndarray
#             Array of length ``n`` containing the canalizing strength of each
#             variable.
#         """
#         if self.n==1:
#             warnings.warn(
#                 "Canalizing strength is only defined for Boolean functions with n > 1 inputs. "
#                 "Returning 1 for n == 1.",
#                 RuntimeWarning
#             )
#             return np.ones(1,dtype=np.float64)
#         res = np.zeros((self.n-1,self.n))
#         for k in range(1, self.n):
#             res[k-1] = self.get_kset_canalizing_proportion_of_variables(k)
#         multipliers = 2**np.arange(1, self.n) / (2**np.arange(1, self.n) - 1)
#         return np.mean(res * multipliers[:, np.newaxis], axis = 0)
        
#     def get_input_redundancy(self) -> float:
#         """
#         Compute the input redundancy of the Boolean function.
    
#         Input redundancy quantifies the fraction of inputs that are not required
#         to determine the output. Constant functions have redundancy 1, whereas
#         parity functions have redundancy 0.
    
#         Returns
#         -------
#         float
#             Normalized input redundancy in the interval ``[0, 1]``.
    
#         Raises
#         ------
#         ImportError
#             If the CANA package is not installed.
    
#         Notes
#         -----
#         This method has exponential time complexity in ``n`` and is intended for
#         small Boolean functions.

#         References
#         ----------
#         Marques-Pita, M., & Rocha, L. M. (2013).
#             Canalization and control in automata networks: body segmentation in
#             *Drosophila melanogaster*. *PLoS One*, 8(3), e55946.
    
#         Correia, R. B., Gates, A. J., Wang, X., & Rocha, L. M. (2018).
#             CANA: a python package for quantifying control and canalization in
#             Boolean networks. *Frontiers in Physiology*, 9, 1046.
#         """
#         utils._require_cana()
#         return self.to_cana().input_redundancy()
        
#     def get_edge_effectiveness(self) -> list[float]:
#         """
#         Compute the edge effectiveness of each input variable.
    
#         Edge effectiveness measures how strongly flipping an input variable
#         influences the output. Non-essential inputs have effectiveness 0, whereas
#         inputs that always flip the output have effectiveness 1.
    
#         Returns
#         -------
#         list[float]
#             List of length ``n`` containing edge effectiveness values in
#             ``[0, 1]``.
    
#         Raises
#         ------
#         ImportError
#             If the CANA package is not installed.

#         Notes
#         -----
#         This method has exponential time complexity in ``n`` and is intended for
#         small Boolean functions.
    
#         References
#         ----------
#         Marques-Pita, M., & Rocha, L. M. (2013).
#             Canalization and control in automata networks: body segmentation in
#             *Drosophila melanogaster*. *PLoS One*, 8(3), e55946.
    
#         Correia, R. B., Gates, A. J., Wang, X., & Rocha, L. M. (2018).
#             CANA: a python package for quantifying control and canalization in
#             Boolean networks. *Frontiers in Physiology*, 9, 1046.
#         """
#         utils._require_cana()
#         return self.to_cana().edge_effectiveness()
    
#     def get_effective_degree(self) -> float:
#         """
#         Compute the effective degree of the Boolean function.
    
#         The effective degree is defined as the sum of the edge effectiveness
#         values of all input variables.
    
#         Returns
#         -------
#         float
#             Effective degree of the Boolean function.
    
#         Raises
#         ------
#         ImportError
#             If the CANA package is not installed.

#         Notes
#         -----
#         This method has exponential time complexity in ``n`` and is intended for
#         small Boolean functions.
        
#         References
#         ----------
#         Marques-Pita, M., & Rocha, L. M. (2013).
#             Canalization and control in automata networks: body segmentation in
#             *Drosophila melanogaster*. *PLoS One*, 8(3), e55946.
    
#         Correia, R. B., Gates, A. J., Wang, X., & Rocha, L. M. (2018).
#             CANA: a python package for quantifying control and canalization in
#             Boolean networks. *Frontiers in Physiology*, 9, 1046.
#         """
#         utils._require_cana()
#         return float(sum(self.get_edge_effectiveness()))



