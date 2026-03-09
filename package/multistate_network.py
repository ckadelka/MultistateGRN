#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import warnings
from collections.abc import Sequence
import numpy as np
import math
from boolforge import WiringDiagram

import utils
from multistate_function import MultistateFunction

class MultistateNetwork(WiringDiagram):
    def __init__(
            self,
            F : Sequence[MultistateFunction | list[int] | np.ndarray],
            R : Sequence[int],
            I : Sequence[Sequence[int]] | WiringDiagram,
            variables : Sequence[str] | None = None
    ):
        if isinstance(F, (str, bytes)) or not isinstance(F, Sequence):
            raise TypeError("F must be a sequence of MultistateFunction objects or truth tables")
        if isinstance(I, (str, bytes)) or not isinstance(I, (Sequence, WiringDiagram)):
            raise TypeError("I must be a sequence of sequences of int or a WiringDiagram instance")
        if isinstance(I, WiringDiagram):
            if variables is not None:
                warnings.warn("Provided variables ignored; using variables from WiringDiagram.", UserWarning)
            super().__init__(I.I, I.variables)
        else:
            super().__init__(I, variables)
        
        if len(F) != self.N:
            raise ValueError("len(F) must match the number of nodes in the wiring diagram")
        if len(F) != len(R):
            raise ValueError("len(F) must match the number of radices len(R)")
        
        self.F = []
        for i, f in enumerate(F):
            if isinstance(f, (Sequence)):
                msf = MultistateFunction(f, R[i], name=self.variables[i])
            elif isinstance(f, MultistateFunction):
                msf = f
                if msf.r != R[i]:
                    raise ValueError("Multistate Function at index {i} has radix {msf.r}, expected {R[i]}")
                msf.name = self.variables[i]
            else:
                raise TypeError(f"Invalid entry in F at index {i}: expected MultistateFunction, truth table, got {type(f)}")
            if msf.n != self.indegrees[i]:
                raise ValueError(f"Index {i}: function has {msf.n} inputs but wiring diagram has degree {self.indegrees[i]}")
            self.F.append(msf)
        
        # duplicate the array
        self.R = []
        for r in R:
            self.R.append(r)
        self.R = np.array(self.R, np.uint8)
        
        self.STG = None

#     def remove_constants(self) -> None:
#         """
#         Remove structurally constant nodes from the Boolean network.
    
#         A node is considered constant if it has no regulators (indegree zero).
#         Such nodes are eliminated from the dynamic network by propagating their
#         fixed Boolean values to downstream nodes. Eliminated constants and their
#         effects are recorded in the ``constants`` attribute.
    
#         Notes
#         -----
#         - The Boolean value of a constant node is taken from its Boolean function.
#         - After removal, ``self.N`` refers to the number of remaining dynamic nodes.
#         - Nodes that lose all regulators as a result of constant removal are
#           assigned a non-essential self-loop to preserve network structure.
#         """
#         # Identify constant nodes from topology
#         # In this model, source nodes (indegree 0) are exactly the semantic constants
#         # at initialization time.
#         indices_constants = self.get_source_nodes(as_dict=False)
#         if len(indices_constants) == 0:
#             return
    
#         dict_constants = self.get_source_nodes(as_dict=True)
#         values_constants = [int(self.F[c][0]) for c in indices_constants]
    
#         # Propagate constant values downstream
#         for id_constant, value in zip(indices_constants, values_constants):    
#             for i in range(self.N):
#                 if dict_constants[i]:
#                     continue
    
#                 try:
#                     index = list(self.I[i]).index(id_constant)
#                 except ValueError:
#                     continue
    
#                 truth_table = utils.get_left_side_of_truth_table(self.indegrees[i])
#                 indices_to_keep = np.where(truth_table[:, index] == value)[0]
    
#                 self.F[i].f = self.F[i].f[indices_to_keep]
    
#                 if self.weights is not None:
#                     self.weights[i] = self.weights[i][self.I[i] != id_constant]
    
#                 self.I[i] = self.I[i][self.I[i] != id_constant]
#                 self.indegrees[i] -= 1
#                 self.F[i].n -= 1
    
    
#             self.constants[str(self.variables[id_constant])] = value
    
#         # Ensure no remaining node loses all regulators
#         for i in range(self.N):
#             if dict_constants[i]:
#                 continue
    
#             if self.indegrees[i] == 0:
#                 self.indegrees[i] = 1
#                 self.F[i].n = 1
#                 self.F[i].f = np.array([self.F[i][0], self.F[i][0]], dtype=int)
#                 self.I[i] = np.array([i], dtype=int)
    
#                 if self.weights is not None:
#                     self.weights[i] = np.array([np.nan], dtype=float)
    
#         # Remove constant nodes structurally (using original mask)
#         self.F = [self.F[i] for i in range(self.N) if not dict_constants[i]]
    
#         adjustment_for_I = np.cumsum([dict_constants[i] for i in range(self.N)])
#         self.I = [
#             self.I[i] - adjustment_for_I[self.I[i]]
#             for i in range(self.N)
#             if not dict_constants[i]
#         ]
    
#         if self.weights is not None:
#             self.weights = [self.weights[i] for i in range(self.N) if not dict_constants[i]]
    
#         self.variables = np.array(
#             [self.variables[i] for i in range(self.N) if not dict_constants[i]],
#             dtype=str,
#         )
    
#         self.indegrees = np.array(
#             [self.indegrees[i] for i in range(self.N) if not dict_constants[i]],
#             dtype=int,
#         )
    
#         # Update network size and recompute outdegrees
#         self.N -= len(indices_constants)
#         self.outdegrees = self.get_outdegrees()

        
#     @classmethod
#     def from_cana(
#         cls,
#         cana_BooleanNetwork: "cana.boolean_network.BooleanNetwork",
#         simplify_functions: bool = False,
#     ) -> "BooleanNetwork":
#         """
#         Construct a BooleanNetwork from a ``cana.BooleanNetwork`` instance.
    
#         This compatibility method converts a Boolean network defined using the
#         ``cana`` package into a BoolForge ``BooleanNetwork``.
    
#         Parameters
#         ----------
#         cana_BooleanNetwork : cana.boolean_network.BooleanNetwork
#             A Boolean network instance from the ``cana`` package.
#         simplify_functions : bool, optional
#             If True, Boolean update functions are simplified after initialization.
#             Default is False.  
            
#         Returns
#         -------
#         BooleanNetwork
#             The corresponding BoolForge BooleanNetwork.
    
#         Raises
#         ------
#         ImportError
#             If the CANA package is not installed.
#         TypeError
#             If the input object does not appear to be a valid CANA BooleanNetwork.
#         KeyError
#             If required fields are missing from the CANA logic specification.
#         """
#         utils._require_cana()
        
#         try:
#             logic = cana_BooleanNetwork.logic
#         except AttributeError as e:
#             raise TypeError(
#                 "Input must be a cana.boolean_network.BooleanNetwork instance."
#             ) from e
    
#         F = []
#         I = []
#         variables = []
    
#         # Ensure deterministic ordering by node index
#         for idx in sorted(logic.keys()):
#             entry = logic[idx]
    
#             if "name" not in entry or "in" not in entry or "out" not in entry:
#                 raise KeyError(
#                     f"Logic entry for node {idx} must contain keys "
#                     "'name', 'in', and 'out'."
#                 )
    
#             variables.append(str(entry["name"]))
#             I.append(list(entry["in"]))
#             F.append(np.array(entry["out"], dtype=int))
    
#         return cls(F=F, I=I, variables=variables, simplify_functions=simplify_functions)

    @classmethod
    def from_string(
        cls,
        network_string: str,
        separator : str = ":",
        original_equality : str = "=",
        original_not : str = "NOT",
        original_and : str = "AND",
        original_or : str = "OR",
        max_degree: int = 24,
        ) -> "MultistateNetwork":
        # """
        # Construct a BooleanNetwork from a textual Boolean rule specification.
    
        # This compatibility method parses a string representation of Boolean update
        # rules (one rule per line) and constructs a corresponding BooleanNetwork. 
        # The input format is intended for legacy or trusted sources and supports 
        # logical expressions using AND/OR/NOT operators. See utils.f_from_expression
        # for details.
    
        # .. warning::
        #     This method uses ``eval`` internally and MUST NOT be used on untrusted
        #     input.
    
        # Parameters
        # ----------
        # network_string : str
        #     String encoding Boolean update rules, one per line.
        # separator : str or sequence of str, optional
        #     Separator(s) between variable names and Boolean expressions.
        # max_degree : int, optional
        #     Maximum allowed indegree for explicit truth-table construction.
        # allow_truncation : bool, optional
        #     If False (default), nodes with indegree greater than ``max_degree``
        #     raise a ValueError. If True, such nodes are replaced by identity
        #     self-loops, allowing fast construction of large networks while
        #     ignoring high-degree functions.
            
        # Returns
        # -------
        # BooleanNetwork
        #     The constructed Boolean network.
    
        # Raises
        # ------
        # ValueError
        #     If parsing fails or if ``allow_truncation`` is False and 
        #     a node exceeds ``max_degree``.
        # """
        
        tvec = network_string
        
        # the first line is inherently unique
        node_count = 1
        
        # get the count of unique nodes -> non-boolean nodes exist across multiple adjacent lines
        for i, line in enumerate(tvec):
            if i > 0 and not line.split(original_equality)[0] == tvec[i - 1].split(original_equality)[0]:
                node_count += 1
        
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
        
        dict_vars_and_consts = dict({original_equality : '==', original_not : '~', original_and : '&', original_or : '|'})
        dict_vars_and_consts.update(dict(list(zip(variables, ["x[%i]" % i for i in range(len(variables))]))))
        # append constants to the end
        dict_vars_and_consts.update(list(zip(constants, ["x[%i]" % i for i in range(len(variables), len(set(constants_and_variables)))])))
        
        # switch tvec to use the new representation instead of the original ones
        for i, line in enumerate(tvec):
            linesplit = line.split(' ')
            for j, el in enumerate(linesplit):
                if el not in ['(', ')', '+', '*', separator, '==', '&', '|', '~', '', ' ']:
                    if original_equality in el:
                        equality_split = el.split(original_equality)
                        linesplit[j] = dict_vars_and_consts[equality_split[0]] + '==' + equality_split[1]
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
                    state_counts = line.split(node_rep + '==')
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
                    state_counts = line.split(node_rep + '==')
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
            regulatee = int(tvec[i][utils.find_all_indices(tvec[i], '[')[0] + 1:utils.find_all_indices(tvec[i], ']')[0]])
            idx_open = utils.find_all_indices(tvec_rule[i], '[')
            idx_close = utils.find_all_indices(tvec_rule[i], ']')
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
            node_idx = int(tvec[i][utils.find_all_indices(tvec[i], '[')[0] + 1:utils.find_all_indices(tvec[i], ']')[0]])
            state = int(line.split(str(node_idx) + ']==')[1].split(' ')[0])
            x = np.zeros(I[node_idx][len(I[node_idx]) - 1] + 1, int)
            for j, _ in enumerate(F[node_idx]):
                vecj = utils.dec2mix(j, B[I[node_idx]])
                for k, I_v in enumerate(I[node_idx]):
                    x[I_v] = vecj[k]
                if eval(tvec_rule[i]):
                    F[node_idx][j] = state
        
        # add self-regulation to F and I for every constant
        for i in range(len(constants)):
            f = []
            for b in range(B[i + len(variables)]):
                f.append(b)
            F.append(np.array(f))
            I.append(np.array([len(variables) + i]))
            degree.append(1)
        
        return cls(F, B, I, variables)

#     @classmethod
#     def from_DiGraph(cls, nx_DiGraph: "nx.DiGraph") -> "WiringDiagram":
#         raise NotImplementedError(
#             "from_DiGraph is not supported for BooleanNetwork. "
#             "Use WiringDiagram.from_DiGraph and then construct "
#             "a BooleanNetwork by providing Boolean update functions."
#         )
    
#     def to_cana(self) -> "cana.boolean_network.BooleanNetwork":
#         """
#         Export the Boolean network as a ``cana.BooleanNetwork`` instance.
    
#         This compatibility method converts the current BooleanNetwork into an
#         equivalent representation from the ``cana`` package. The exported network
#         reflects the current state of the model, including any removed constants,
#         simplifications, or identity self-loops.
    
#         Returns
#         -------
#         cana.boolean_network.BooleanNetwork
#             A ``cana`` BooleanNetwork instance representing this network.
    
#         Raises
#         ------
#         ImportError
#             If the ``cana`` package is not installed.
#         """
#         try:
#             import cana.boolean_network
#         except ImportError as e:
#             raise ImportError(
#                 "The 'cana' package is required for to_cana()."
#             ) from e
    
#         logic_dicts = []
#         for bf, regulators, var in zip(self.F, self.I, self.variables):
#             logic_dicts.append(
#                 {
#                     "name": var,
#                     "in": list(regulators),
#                     "out": bf.f.tolist(),
#                 }
#             )
    
#         return cana.boolean_network.BooleanNetwork(
#             Nnodes=self.N,
#             logic={i: d for i, d in enumerate(logic_dicts)},
#         )


#     def to_bnet(
#         self,
#         separator: str = ",\t",
#         as_polynomial: bool = True,
#     ) -> str:
#         """
#         Export the Boolean network in BNET format.
    
#         This compatibility method returns a string representation of the Boolean
#         network in the BNET format used by tools such as BoolNet and PyBoolNet,
#         with one line per variable of the form ``variable <separator> function.
        
#         Parameters
#         ----------
#         separator : str, optional
#             String used to separate the target variable from its update function.
#             Default is `",\\t"`.
#         as_polynomial : bool, optional
#             If True (default), return Boolean functions in polynomial form.
#             If False, return functions as logical expressions.
    
#         Returns
#         -------
#         str
#             A string containing the BNET representation of the network.
            
#         Notes
#         -----
#         This method exports the reduced Boolean network, i.e. after semantic
#         constants have been removed during initialization.
#         """
#         lines = []
    
#         for i in range(self.N):
#             if as_polynomial:
#                 function = utils.bool_to_poly(
#                     self.F[i],
#                     self.variables[self.I[i]].tolist(),
#                 )
#             else:
#                 function = self.F[i].to_expression(" & ", " | ")
    
#             lines.append(f"{self.variables[i]}{separator}{function}")
    
#         return "\n".join(lines)
    
    
#     def to_truth_table(
#         self,
#         filename: str | None = None,
#     ) -> pd.DataFrame:
#         """
#         Construct the full synchronous truth table of the Boolean network.
    
#         Each row corresponds to a network state at time ``t`` and its deterministic
#         successor at time ``t+1`` under synchronous updating.
    
#         Parameters
#         ----------
#         filename : str, optional
#             If provided, the truth table is written to a file. The file extension
#             determines the format and must be one of ``'csv'``, ``'xls'``, or
#             ``'xlsx'``. If None (default), no file is created.
    
#         Returns
#         -------
#         pandas.DataFrame
#             The full truth table with shape ``(2**N, 2*N)``.
    
#         Notes
#         -----
#         - States are enumerated in lexicographic order, consistent with
#           ``utils.get_left_side_of_truth_table``.
#         - This method computes and stores the synchronous state transition graph
#           (``self.STG``) if it has not been computed previously.
#         - Exporting to Excel requires the ``openpyxl`` package.
#         """
#         columns = [name + "(t)" for name in self.variables]
#         columns += [name + "(t+1)" for name in self.variables]
    
#         if self.STG is None:
#             self.compute_synchronous_state_transition_graph()
    
#         data = np.zeros((2**self.N, 2*self.N), dtype=int)
#         data[:, :self.N] = utils.get_left_side_of_truth_table(self.N)
    
#         for i in range(2**self.N):
#             data[i, self.N:] = utils.dec2bin(self.STG[i], self.N)
    
#         truth_table = pd.DataFrame(data, columns=columns)
    
#         if filename is not None:
#             if not isinstance(filename, str):
#                 raise TypeError("filename must be a string")
    
#             ending = filename.split(".")[-1]
#             if ending not in {"csv", "xls", "xlsx"}:
#                 raise ValueError("filename must end in 'csv', 'xls', or 'xlsx'")
    
#             if ending == "csv":
#                 truth_table.to_csv(filename, index=False)
#             else:
#                 truth_table.to_excel(filename, index=False)
    
#         return truth_table

    def __len__(self):
        return self.N
    
    def __str__(self):
        return f"MultistateNetwork(N={self.N}, indegrees={self.indegrees.tolist()})"
    
    def __getitem__(self, index):
        return self.F[index]
    
    def __repr__(self):
        return f"{type(self).__name__}(N={self.N})"
    
    def __call__(self, state):
        """
        Apply one synchronous update step to the Multstate network.
        
        The next state is obtained by evaluating each node's multistate update
        function on the current values of its regulators.

        Parameters
        ----------
        state : Sequence of int
            Current network state as a mixed-radix vector of length ``N``, ordered
            according to ``self.variables``.

        Returns
        -------
        np.ndarray
            The updated network state after one synchronous update.
        
        Notes
        -----
        This method is equivalnet to calling ``update_network_synchronosly``.
        """
        return self.update_network_synchronously(state)
    
    
#     def get_types_of_regulation(self) -> list[np.ndarray]:
#         """
#         Compute and return regulation types (weights) for all nodes in the network.
    
#         For each Boolean function, the type of each input regulation is determined
#         via ``BooleanFunction.get_type_of_inputs`` and mapped to numerical weights
#         using ``dict_weights``. The resulting weights are stored in the
#         ``self.weights`` attribute and also returned.
    
#         Returns
#         -------
#         list of np.ndarray
#             Regulation weights for each node, aligned with the wiring diagram.
    
#         Notes
#         -----
#         - This method recomputes ``self.weights`` from scratch.
#         - Calling this method overwrites any existing values in ``self.weights``.
#         """
#         self.weights = [
#             np.array([dict_weights[el] for el in bf.get_type_of_inputs()], dtype=float)
#             for bf in self.F
#         ]
#         return self.weights



#     ## Transform Boolean networks
#     def simplify_functions(self) -> None:
#         """
#         Remove all non-essential regulators from the Boolean network.
    
#         For each node, non-essential regulators (identified via ``np.nan`` entries
#         in ``self.weights``) are removed from the wiring diagram and the associated
#         Boolean function is restricted to its essential inputs. Nodes that would
#         otherwise lose all regulators are assigned an identity self-loop to preserve
#         network structure.
    
#         Notes
#         -----
#         - This method modifies the network in place.
#         - Regulation types (``self.weights``) are recomputed if necessary.
#         - Identity self-loops introduced here are structural artifacts and do not
#           represent genuine regulatory interactions.
#         """
#         # Ensure regulation types / weights are available
#         self.get_types_of_regulation()
    
#         for i in range(self.N):
#             regulator_is_non_essential = np.isnan(self.weights[i])
    
#             # All regulators are essential
#             if not np.any(regulator_is_non_essential):
#                 continue
    
#             non_essential_variables = np.where(regulator_is_non_essential)[0]
#             essential_variables = np.where(~regulator_is_non_essential)[0]
    
#             # Update outdegrees (each regulator appears at most once in I[i])
#             self.outdegrees[non_essential_variables] -= 1
    
#             # No essential regulators: introduce identity self-loop
#             if len(essential_variables) == 0:
#                 self.indegrees[i] = 1
#                 self.F[i].f = np.array([self.F[i][0], self.F[i][0]], dtype=int)
#                 self.F[i].n = 1
#                 self.F[i].variables = np.array([self.variables[i]], dtype=str)
#                 self.I[i] = np.array([i], dtype=int)
#                 self.weights[i] = np.array([np.nan], dtype=float)
#                 self.outdegrees[i] += 1  # keep sum(outdegrees) == sum(indegrees)
#                 continue
    
#             # Restrict truth table to essential inputs
#             left_side = utils.get_left_side_of_truth_table(self.indegrees[i])
#             mask = np.sum(left_side[:, non_essential_variables], axis=1) == 0
    
#             self.F[i].f = self.F[i][mask]
#             self.F[i].n = len(essential_variables)
#             self.F[i].variables = self.F[i].variables[~regulator_is_non_essential]
#             self.I[i] = self.I[i][essential_variables]
#             self.weights[i] = self.weights[i][essential_variables]
#             self.indegrees[i] = len(essential_variables)


#     def get_identity_nodes(
#         self, 
#         as_dict: bool = False
#     ) -> dict[int, bool] | np.ndarray:
#         """
#         Identify identity (memory) nodes in the Boolean network.
    
#         An identity node is a node with a single self-regulatory edge whose
#         Boolean update function is the identity function ``f(x) = x``. Such
#         nodes retain their state over time unless externally modified.
    
#         Parameters
#         ----------
#         as_dict : bool, optional
#             If True, return a dictionary mapping node indices to booleans.
#             If False (default), return an array of indices of identity nodes.
    
#         Returns
#         -------
#         dict[int, bool] or np.ndarray
#             If ``as_dict`` is True, a dictionary indicating which nodes are
#             identity nodes.
#             If ``as_dict`` is False, an array of indices of identity nodes.
#         """
#         is_identity = np.array(
#             [
#                 self.indegrees[i] == 1
#                 and self.I[i][0] == i
#                 and self.F[i][0] == 0
#                 and self.F[i][1] == 1
#                 for i in range(self.N)
#             ],
#             dtype=bool,
#         )
    
#         if as_dict:
#             return dict(enumerate(is_identity.tolist()))
#         return np.where(is_identity)[0]
    
    
#     def propagate_constants(self) -> "BooleanNetwork":
#         """
#         Recursively propagate constants through the network.
    
#         Any node whose update function becomes constant is converted
#         into a structural constant. Removal of such nodes and updating
#         of self.constants is handled by __init__.
#         """
    
#         F = deepcopy(self.F)
#         I = deepcopy(self.I)
    
#         n = len(F)
    
#         # Build reverse dependency graph
#         dependents = {i: [] for i in range(n)}
#         for node in range(n):
#             for inp in I[node]:
#                 dependents[inp].append(node)
    
#         fixed = {}
#         queue = deque()
#         indices_fixation_layers = []  # <-- new
    
#         # ----------------------------------------------------------
#         # Initial scan (Layer 0)
#         # ----------------------------------------------------------
#         initial_layer = []
    
#         for node in range(n):
#             if F[node].is_constant():
#                 constant_value = F[node][0]
#                 fixed[node] = constant_value
#                 queue.append(node)
#                 initial_layer.append(node)
    
#                 F[node].f = np.array([constant_value], dtype=int)
#                 F[node].n = 0
#                 I[node] = np.array([], dtype=int)
    
#         if initial_layer:
#             indices_fixation_layers.append(initial_layer)
    
#         # ----------------------------------------------------------
#         # Propagation
#         # ----------------------------------------------------------
#         while queue:
    
#             current_layer_size = len(queue)
#             next_layer = []
    
#             for _ in range(current_layer_size):
    
#                 fixed_node = queue.popleft()
#                 value = fixed[fixed_node]
    
#                 for node in dependents[fixed_node]:
    
#                     if node in fixed:
#                         continue
    
#                     inputs = I[node]
#                     positions = np.where(inputs == fixed_node)[0]
#                     if len(positions) == 0:
#                         continue
    
#                     pos = positions[0]
#                     k = len(inputs)
#                     old_table = F[node].f
#                     new_table = []
    
#                     for r in range(len(old_table)):
#                         bit = (r >> (k - pos - 1)) & 1
#                         if bit == value:
#                             new_table.append(old_table[r])
    
#                     new_table = np.array(new_table, dtype=int)
    
#                     F[node].f = new_table
#                     F[node].n = k - 1
#                     I[node] = np.delete(inputs, pos)
    
#                     if len(new_table) == 1 or np.all(new_table == new_table[0]):
#                         const_value = int(new_table[0])
    
#                         F[node].f = np.array([const_value], dtype=int)
#                         F[node].n = 0
#                         I[node] = np.array([], dtype=int)
    
#                         fixed[node] = const_value
#                         queue.append(node)
#                         next_layer.append(node)
    
#             if next_layer:
#                 indices_fixation_layers.append(next_layer)
        
#         fixation_layers = [
#             [str(self.variables[i]) for i in layer]
#             for layer in indices_fixation_layers
#         ]
        
#         # Reinitialize — constructor removes structural constants
#         reduced_bn = self.__class__(F, I, self.variables)
#         return {'ReducedNetwork' : reduced_bn,
#                 'FixationLayers' : fixation_layers}


#     def get_network_with_fixed_identity_nodes(
#         self,
#         values_identity_nodes: Sequence[int],
#         keep_controlled_nodes: bool = False,
#         simplify_recursively: bool = False,
#     ) -> "BooleanNetwork":
#         """
#         Construct a Boolean network with identity nodes fixed to given values.
    
#         Identity nodes are nodes with a single self-regulatory edge and identity
#         update rule ``f(x) = x``. This method fixes the values of such nodes and
#         returns a new BooleanNetwork with the corresponding constants removed.
    
#         Parameters
#         ----------
#         values_identity_nodes : sequence of int
#             Values to fix for each identity node, in the order returned by
#             ``get_identity_nodes(as_dict=False)``. Each value must be either 0 or 1.
#         keep_controlled_nodes : bool, optional
#             If True, controlled nodes are retained in the network as identity
#             nodes with self-loops. If False (default), controlled nodes are
#             eliminated as constants.
#         simplify_recursively : bool, optional
#             If True, recursively propagate fixed values through the network and
#             eliminate any nodes whose update functions become constant as a result.
#             This logical reduction is repeated until no further simplifications are
#             possible. If False (default), only the explicitly controlled nodes are
#             fixed and no additional recursive simplification is performed. 
            
#         Returns
#         -------
#         BooleanNetwork
#             A new BooleanNetwork with the specified identity nodes fixed.
#         """
#         indices_identity_nodes = self.get_identity_nodes(as_dict=False)
    
#         if len(values_identity_nodes) != len(indices_identity_nodes):
#             raise ValueError(
#                 f"The number of values provided ({len(values_identity_nodes)}) must "
#                 f"match the number of identity nodes ({len(indices_identity_nodes)})."
#             )
    
#         for v in values_identity_nodes:
#             if v not in (0, 1):
#                 raise ValueError("Identity node values must be 0 or 1.")
    
#         return self.get_network_with_node_controls(
#             indices_controlled_nodes=indices_identity_nodes,
#             values_controlled_nodes=values_identity_nodes,
#             keep_controlled_nodes=keep_controlled_nodes,
#             simplify_recursively=simplify_recursively,
#         )


#     def get_network_with_node_controls(
#         self,
#         indices_controlled_nodes: Sequence[int],
#         values_controlled_nodes: Sequence[int],
#         keep_controlled_nodes: bool = False,
#         simplify_recursively: bool = False,
#     ) -> "BooleanNetwork":
#         """
#         Construct a Boolean network with specified nodes fixed to given values.
    
#         This method applies node-level interventions by fixing selected nodes to
#         constant Boolean values. Controlled nodes may either be removed from the
#         dynamic network as constants or retained as identity-clamped nodes.
    
#         Parameters
#         ----------
#         indices_controlled_nodes : sequence of int
#             Indices of nodes to be fixed.
#         values_controlled_nodes : sequence of int
#             Values to fix for each specified node, in the same order as
#             ``indices_controlled_nodes``. Each value must be either 0 or 1.
#         keep_controlled_nodes : bool, optional
#             If True, controlled nodes are retained in the network as identity
#             nodes with self-loops. If False (default), controlled nodes are
#             eliminated as constants.
#         simplify_recursively : bool, optional
#             If True, recursively propagate fixed values through the network and
#             eliminate any nodes whose update functions become constant as a result.
#             This logical reduction is repeated until no further simplifications are
#             possible. If False (default), only the explicitly controlled nodes are
#             fixed and no additional recursive simplification is performed. 
            
#         Returns
#         -------
#         BooleanNetwork
#             A new BooleanNetwork with the specified node controls applied.
#         """
#         if len(indices_controlled_nodes) != len(values_controlled_nodes):
#             raise ValueError(
#                 f"The number of controlled nodes ({len(indices_controlled_nodes)}) "
#                 f"must match the number of values provided ({len(values_controlled_nodes)})."
#             )
    
#         for node in indices_controlled_nodes:
#             if not isinstance(node, (int, np.integer)) or node < 0 or node >= self.N:
#                 raise ValueError(f"Invalid node index: {node, not isinstance(node, int), node < 0 , node >= self.N}")
    
#         for v in values_controlled_nodes:
#             if v not in (0, 1):
#                 raise ValueError("Controlled node values must be 0 or 1.")
    
#         if simplify_recursively and keep_controlled_nodes:
#             raise ValueError(
#                 "Cannot simplify recursively when keep_controlled_nodes=True."
#             )
    
#         F = deepcopy(self.F)
#         I = deepcopy(self.I)
#         controlled_variables = [str(self.variables[int(i)]) for i in indices_controlled_nodes]
    
#         for node, value in zip(indices_controlled_nodes, values_controlled_nodes):
#             if keep_controlled_nodes:
#                 # Identity-clamped node
#                 F[node].f = np.array([value, value], dtype=int)
#                 I[node] = np.array([node], dtype=int)
#             else:
#                 # Structural constant (to be removed)
#                 F[node].f = np.array([value], dtype=int)
#                 F[node].n = 0
#                 I[node] = np.array([], dtype=int)
    
#         bn2 = self.__class__(F, I, self.variables) #__init__ removes fixated control nodes
    
#         # Preserve previously removed constants if controlled nodes are eliminated
#         if simplify_recursively:
#             dummy = bn2.propagate_constants()
#             bn3 = dummy['ReducedNetwork']
#             fixation_layers = dummy['FixationLayers']
#             bn3.constants.update(self.constants)
#             bn3.constants.update(bn2.constants)
#             bn3.fixation_layers = [controlled_variables] + fixation_layers
#             return bn3
#         else:
#             if not keep_controlled_nodes:
#                 bn2.constants.update(self.constants)
#                 bn2.fixation_layers = [controlled_variables]
#             return bn2

#     def get_network_with_edge_controls(
#         self,
#         control_targets: Sequence[int],
#         control_sources: Sequence[int],
#         values_edge_controls: Sequence[int] | None = None,
#         keep_fully_controlled_nodes: bool = True
#     ) -> "BooleanNetwork":
#         """
#         Construct a Boolean network with specified regulatory edges controlled.
    
#         This method fixes the influence of selected source nodes on selected target
#         nodes by restricting the target's Boolean update function to entries where
#         the source assumes a specified value, and then removing the corresponding
#         regulatory edge.
    
#         Parameters
#         ----------
#         control_targets : sequence of int
#             Indices of target nodes.
#         control_sources : sequence of int
#             Indices of source nodes whose influence on the corresponding targets
#             is to be controlled.
#         values_edge_controls : sequence of int, optional
#             Fixed values (0 or 1) imposed on each controlled edge. If None, all
#             controlled edges are fixed to 0.
#         keep_fully_controlled_nodes : bool, optional
#             If True (default), nodes without any remaining regulation are retained
#             in the network as identity nodes with self-loops. 
#             If False, fully controlled nodes are eliminated as constants.
            
#         Returns
#         -------
#         BooleanNetwork
#             A new BooleanNetwork with the specified edge controls applied.
    
#         Raises
#         ------
#         ValueError
#             If input lengths do not match, indices are invalid, or edge values are
#             not in {0, 1}.
#         """
#         if len(control_targets) != len(control_sources):
#             raise ValueError("control_targets and control_sources must have equal length.")
    
#         if values_edge_controls is None:
#             values_edge_controls = [0] * len(control_targets)
    
#         if len(values_edge_controls) != len(control_targets):
#             raise ValueError(
#                 "values_edge_controls must have the same length as control_targets."
#             )
    
#         F_new = deepcopy(self.F)
#         I_new = deepcopy(self.I)
    
#         for target, source, fixed_value in zip(
#             control_targets, control_sources, values_edge_controls
#         ):
#             if fixed_value not in (0, 1):
#                 raise ValueError("Edge control values must be 0 or 1.")
    
#             if not (0 <= target < self.N):
#                 raise ValueError(f"Invalid target index: {target}")
    
#             if not (0 <= source < self.N):
#                 raise ValueError(f"Invalid source index: {source}")
    
#             if source not in I_new[target]:
#                 raise ValueError(
#                     f"Source node {source} is not a regulator of target node {target}."
#                 )
    
#             idx_reg = list(I_new[target]).index(source)
#             n_inputs = F_new[target].n
    
#             truth_indices = np.arange(2**n_inputs, dtype=np.uint32)
#             mask = ((truth_indices >> (n_inputs - 1 - idx_reg)) & 1) == fixed_value
    
#             # Restrict truth table
#             F_new[target].f = F_new[target].f[mask]
#             F_new[target].n -= 1
    
#             # Remove regulator
#             I_new[target] = np.delete(I_new[target], idx_reg)
            
#             # ---- NEW LOGIC: fully controlled node -----------------------------
#             if F_new[target].n == 0:
#                 if keep_fully_controlled_nodes:
#                     value = int(F_new[target].f[0])
#                     # Identity-clamped node
#                     F_new[target].f = np.array([value, value], dtype=int)
#                     F_new[target].n = 1
#                     I_new[target] = np.array([target], dtype=int)

    
#         return self.__class__(F_new, I_new, self.variables)

#def update_ms_single_node(F, states_regulators, B):
#    return F[multi2dec(states_regulators, B)]
        
    def update_single_node(self, index : int, regulators : Sequence[int]) -> int:
        """
        Update the state of a single node.

        Parameters
        ----------
        index : int
            Index of the node to update.
        regulators : Sequence[int]
            Mixed-radix states of the node's regulators.
        radices : Sequence[int]
            The radices of the provided regulators.

        Returns
        -------
        int
            Updated state of the node.
        """
        return self.F[index].f[utils.mix2dec(regulators, self.F[index].in_rs)].item()
    
    def update_network_synchronously(self, state : Sequence[int]) -> np.ndarray:
        """
        Perform a synchronus update of the multistate network.

        Parameters
        ----------
        state : Sequence[int]
            Mixed-radix state vector of length ``N``.

        Returns
        -------
        np.ndarray
            Updated state vector.

        """
        state = np.asarray(state, int)
        if state.shape[0] != self.N:
            raise ValueError(f"State vector must have length {self.N}, got {state.shape[0]}.")
        if not np.all(state < self.R):
            raise ValueError("State vector cannot exceed defined radices.")
        return self._update_network_synchronously_unchecked(state)
    
    def _update_network_synchronously_unchecked(self, state : np.ndarray) -> np.ndarray:
        next_state = np.zeros(self.N, int)
        for i in range(self.N):
            next_state[i] = self.update_ms_single_node(state[self.I[i]])
        return next_state
    
#     def update_network_SDDS(
#         self,
#         state: Sequence[int],
#         P: np.ndarray,
#         *,
#         rng=None,
#     ) -> np.ndarray:
#         """
#         Perform a stochastic discrete dynamical system (SDDS) update of the network.
    
#         This update scheme follows the SDDS formalism: for each node, the
#         deterministic Boolean update is first computed. If the update would
#         increase the node's state, the change occurs with the node-specific
#         activation probability. If the update would decrease the node's state,
#         the change occurs with the node-specific degradation probability.
#         Otherwise, the node's state remains unchanged.
    
#         Parameters
#         ----------
#         state : sequence of int
#             Current network state (binary vector of length ``N``).
#         P : np.ndarray
#             Array of shape ``(N, 2)``, where ``P[i, 0]`` is the activation
#             probability and ``P[i, 1]`` is the degradation probability for node ``i``.
#         rng : optional
#             Random number generator or seed, passed to ``utils._coerce_rng``.
    
#         Returns
#         -------
#         np.ndarray
#             Updated network state after one stochastic SDDS update.
    
#         Notes
#         -----
#         This implementation follows the SDDS framework introduced in:
    
#         Murrugarra, D., Veliz-Cuba, A., Aguilar, B., Arat, S., & Laubenbacher, R.
#         (2012). *Modeling stochasticity and variability in gene regulatory networks*.
#         EURASIP Journal on Bioinformatics and Systems Biology, 2012(1), 5.
    
#         The method assumes that ``state`` is a valid binary vector and that
#         ``P`` has the correct shape; no additional validation is performed
#         for performance reasons.
#         """
#         rng = utils._coerce_rng(rng)
#         state = np.asarray(state, dtype=int)
    
#         Fx = state.copy()
#         for i in range(self.N):
#             nextstep = self.update_single_node(
#                 index=i,
#                 states_regulators=state[self.I[i]],
#             )
    
#             if nextstep > state[i]:
#                 if rng.random() < P[i, 0]:
#                     Fx[i] = nextstep
#             elif nextstep < state[i]:
#                 if rng.random() < P[i, 1]:
#                     Fx[i] = nextstep
    
#         return Fx



#     def get_steady_states_asynchronous_exact(
#         self,
#         stochastic_weights: Sequence[float] | None = None,
#         max_iterations: int = 1000,
#         tol: float = 1e-9,
#     ) -> dict:
#         """
#         Compute the steady states and basin probabilities under general asynchronous update.
        
#         This method exhaustively constructs the asynchronous state transition graph
#         (STG) of the Boolean network under a general asynchronous update scheme,
#         where nodes are selected for update according to given propensities.
#         The resulting Markov chain is solved exactly using an iterative
#         Gauss–Seidel scheme to obtain absorption probabilities into steady states.
        
#         Parameters
#         ----------
#         stochastic_weights : sequence of float or None, optional
#             Relative update propensities for each node. If None (default),
#             all nodes are updated with equal probability. The weights are
#             normalized internally.
#         max_iterations : int, optional
#             Maximum number of Gauss–Seidel iterations used to compute absorption
#             probabilities before declaring non-convergence.
#         tol : float, optional
#             Convergence tolerance for the infinity norm of probability updates.
#         s
        
#         Returns
#         -------
#         dict
#             Dictionary with the following entries:
        
#             - ``"SteadyStates"`` : list of int  
#               Decimal representations of steady states.
        
#             - ``"NumberOfSteadyStates"`` : int  
#               Total number of steady states.
        
#             - ``"BasinSizes"`` : np.ndarray  
#               Fraction of the state space converging to each steady state.
        
#             - ``"STGAsynchronous"`` : dict  
#               Asynchronous state transition graph represented as a Markov kernel.
        
#             - ``"FinalTransitionProbabilities"`` : np.ndarray  
#               Absorption probabilities from each state into each steady state.
        
#         Raises
#         ------
#         RuntimeError
#             If the iterative solver does not converge within ``max_iterations``.
#         """

#         left_side_of_truth_table = utils.get_left_side_of_truth_table(self.N)

#         if stochastic_weights is not None:
#             stochastic_weights = np.asarray(stochastic_weights, dtype=float)
#             if stochastic_weights.shape[0] != self.N:
#                 raise ValueError("stochastic_weights must have length N.")
#             if np.any(stochastic_weights <= 0):
#                 raise ValueError("stochastic_weights must be strictly positive.")
#             stochastic_weights = stochastic_weights / stochastic_weights.sum()
#         else:
#             stochastic_weights = np.full(self.N, 1.0 / self.N)
            
#         steady_states = []
#         steady_state_dict = {}
#         STG = dict(zip(range(2**self.N),[{} for i in range(2**self.N)]))
#         sped_up_STG = dict(zip(range(2**self.N),[[np.zeros(0,dtype=int),np.zeros(0,dtype=float)] for i in range(2**self.N)]))
#         for xdec in range(2**self.N):
#             x = left_side_of_truth_table[xdec].copy() #important: must create a copy here!
#             to_be_distributed = 0
#             for i in range(self.N):
#                 fx_i = self.update_single_node(i, x[self.I[i]])
#                 if fx_i > x[i]:
#                     fxdec = xdec + 2**(self.N - 1 - i)
#                 elif fx_i < x[i]:
#                     fxdec = xdec - 2**(self.N - 1 - i)
#                 else:
#                     fxdec = xdec
#                 if fxdec in STG[xdec]:
#                     STG[xdec][fxdec] += float(stochastic_weights[i])
#                 else:
#                     STG[xdec][fxdec] = float(stochastic_weights[i])
#                 if fxdec!=xdec:
#                     sped_up_STG[xdec][0] = np.append(sped_up_STG[xdec][0], fxdec)
#                     sped_up_STG[xdec][1] = np.append(sped_up_STG[xdec][1], stochastic_weights[i])
#                 else:
#                     to_be_distributed += stochastic_weights[i]
#             if to_be_distributed < 1:
#                 sped_up_STG[xdec][1] /= (1-to_be_distributed)
#             if len(STG[xdec])==1:
#                 steady_state_dict[xdec] = len(steady_states)
#                 steady_states.append(xdec)
#                 sped_up_STG[xdec][0] = np.append(sped_up_STG[xdec][0], xdec)
#                 sped_up_STG[xdec][1] = np.append(sped_up_STG[xdec][1], 1)
                
#         # Probability vectors for all states
#         final_probabilities = np.zeros((2**self.N, len(steady_states)), dtype=float)
    
#         # Boundary conditions: absorbing states have probability 1 of themselves
#         for xdec in steady_states:
#             final_probabilities[xdec, steady_state_dict[xdec]] = 1.0
#         transient_states = [xdec for xdec in range(2**self.N) if xdec not in steady_state_dict]
        
#         for it in range(1, max_iterations + 1):
#             max_delta = 0.0
    
#             # In-place Gauss–Seidel  update:
#             for xdec in transient_states:
#                 nxt, pr = sped_up_STG[xdec]

#                 old = final_probabilities[xdec].copy()
#                 final_probabilities[xdec] = np.dot(pr, final_probabilities[nxt, :])   # weighted average of successor probability vectors
    
#                 # track convergence (infinity norm per row)
#                 delta = np.max(np.abs(final_probabilities[xdec] - old))
#                 if delta > max_delta:
#                     max_delta = delta
    
#             if max_delta < tol:
#                 basin_sizes = final_probabilities.sum(0)/2**self.N
                
#                 return {
#                     "SteadyStates": steady_states,
#                     "NumberOfSteadyStates": len(steady_states),
#                     "BasinSizes": basin_sizes,
#                     "STGAsynchronous": STG,
#                     "FinalTransitionProbabilities": final_probabilities,
#                 }
            
#         raise RuntimeError(f"Did not converge in {max_iterations} iterations; last max_delta={max_delta:g}")
        

#     def get_steady_states_asynchronous(
#         self,
#         n_simulations: int = 500,
#         initial_states: Sequence[int] | None = None,
#         search_depth: int = 50,
#         debug: bool = False,
#         *,
#         rng=None,
#     ) -> dict:
#         """
#         Approximate steady states of a Boolean network under asynchronous updates.
    
#         This method performs a Monte Carlo–style exploration of the asynchronous
#         state space by simulating asynchronous updates from a collection of initial
#         states. Each simulation proceeds until a steady state is reached or until
#         a maximum search depth is exceeded.
    
#         Unlike ``get_steady_states_asynchronous_exact``, this method does *not*
#         exhaustively explore the full state space and does not guarantee that all
#         steady states will be found. It is intended for large networks where exact
#         enumeration is infeasible.
    
#         Parameters
#         ----------
#         n_simulations : int, optional
#             Number of asynchronous simulations to perform (default is 500).
#         initial_states : sequence of int or None, optional
#             Initial states to use for the simulations, given as decimal
#             representations of network states. If None (default), ``n_simulations``
#             random initial states are generated.
#         search_depth : int, optional
#             Maximum number of asynchronous update steps per simulation before
#             giving up on convergence (default is 50).
#         debug : bool, optional
#             If True, print detailed debugging information during simulation.
#         rng : optional
#             Random number generator or seed, passed to ``utils._coerce_rng``.
    
#         Returns
#         -------
#         dict
#             Dictionary with the following entries:
    
#             - ``"SteadyStates"`` : list of int  
#               Decimal representations of steady states encountered.
    
#             - ``"NumberOfSteadyStates"`` : int  
#               Number of unique steady states found.
    
#             - ``"BasinSizes"`` : list of int  
#               Fraction of simulations converged to each steady state.
    
#             - ``"STGAsynchronous"`` : dict  
#               Partial cache of asynchronous transitions encountered during
#               simulation. Keys are ``(state, node_index)`` and values are
#               successor states (all in decimal form).
    
#             - ``"InitialSamplePoints"`` : list of int  
#               Decimal initial states used in the simulations (either provided
#               explicitly or generated randomly).
    
#         Notes
#         -----
#         - This method detects only *steady states* (fixed points). If the
#           asynchronous dynamics contain limit cycles, simulations may fail
#           to converge within ``search_depth``.
#         - The returned asynchronous transition graph is generally incomplete
#           and should be interpreted as a cache of explored transitions rather
#           than the full STG.
#         """
#         rng = utils._coerce_rng(rng)
    
#         sampled_states: list[int] = []
#         STG_asynchronous: dict[tuple[int, int], int] = {}
    
#         steady_states: list[int] = []
#         basin_sizes: list[int] = []
#         steady_state_dict: dict[int, int] = {}
    
#         for iteration in range(n_simulations):
#             # Initialize state
#             if initial_states is None:
#                 x = rng.integers(2, size=self.N)
#                 xdec = utils.bin2dec(x)
#                 sampled_states.append(xdec)
#             else:
#                 xdec = initial_states[iteration]
#                 x = utils.dec2bin(xdec, self.N)
    
#             if debug:
#                 print(iteration, -1, -1, False, xdec, x)
    
#             for step in range(search_depth):
#                 found_new_state = False
    
#                 # Check if state is already known to be steady
#                 if xdec in steady_state_dict:
#                     basin_sizes[steady_state_dict[xdec]] += 1
#                     break
    
#                 update_order = rng.permutation(self.N)
#                 for i in map(int, update_order):
#                     try:
#                         fxdec = STG_asynchronous[(xdec, i)]
#                     except KeyError:
#                         fx_i = self.update_single_node(i, x[self.I[i]])
#                         if fx_i > x[i]:
#                             fxdec = xdec + 2 ** (self.N - 1 - i)
#                             x[i] = 1
#                             found_new_state = True
#                         elif fx_i < x[i]:
#                             fxdec = xdec - 2 ** (self.N - 1 - i)
#                             x[i] = 0
#                             found_new_state = True
#                         else:
#                             fxdec = xdec
#                         STG_asynchronous[(xdec, i)] = fxdec
    
#                     if fxdec != xdec:
#                         xdec = fxdec
#                         found_new_state = True
#                         break
    
#                 if debug:
#                     print(iteration, step, i, found_new_state, xdec, x)
    
#                 if not found_new_state:
#                     # New steady state found
#                     if xdec in steady_state_dict:
#                         basin_sizes[steady_state_dict[xdec]] += 1
#                     else:
#                         steady_state_dict[xdec] = len(steady_states)
#                         steady_states.append(xdec)
#                         basin_sizes.append(1)
#                     break
    
#             if debug:
#                 print()
    
#         if sum(basin_sizes) < n_simulations:
#             print(
#                 f"Warning: only {sum(basin_sizes)} of the {n_simulations} simulations "
#                 "reached a steady state. Consider increasing search_depth. "
#                 "The network may also contain asynchronous limit cycles."
#             )
        
#         if sum(basin_sizes)>0:
#             basin_sizes /= sum(basin_sizes)
        
#         return {
#             "SteadyStates": steady_states,
#             "NumberOfSteadyStates": len(steady_states),
#             "BasinSizes": basin_sizes,
#             "STGAsynchronous": STG_asynchronous,
#             "InitialSamplePoints": (
#                 initial_states if initial_states is not None else sampled_states
#             ),
#         }


#     def get_steady_states_asynchronous_given_one_initial_condition(
#         self,
#         initial_condition: int | Sequence[int] = 0,
#         n_simulations: int = 500,
#         stochastic_weights: Sequence[float] | None = None,
#         search_depth: int = 50,
#         debug: bool = False,
#         *,
#         rng=None,
#     ) -> dict:
#         """
#         Approximate steady states reachable from a single initial condition under
#         asynchronous updates.
    
#         This method performs multiple asynchronous simulations starting from the
#         same initial condition. In each simulation, nodes are updated one at a time
#         according to either a uniform random order or node-specific stochastic
#         update propensities. The simulation proceeds until a steady state is reached
#         or a maximum number of update steps is exceeded.
    
#         The method is sampling-based and does *not* guarantee that all reachable
#         steady states are found. It is intended for exploratory analysis and for
#         networks where exhaustive asynchronous analysis is infeasible.
    
#         Parameters
#         ----------
#         initial_condition : int or sequence of int, optional
#             Initial network state. If an integer is provided, it is interpreted as
#             the decimal encoding of a Boolean state. If a sequence is provided, it
#             must be a binary vector of length ``N``. Default is 0.
#         n_simulations : int, optional
#             Number of asynchronous simulation runs (default is 500).
#         stochastic_weights : sequence of float or None, optional
#             Relative update propensities for each node. If provided, must have
#             length ``N`` and be strictly positive. The weights are normalized
#             internally. If None (default), nodes are updated uniformly at random.
#         search_depth : int, optional
#             Maximum number of asynchronous update steps per simulation.
#         debug : bool, optional
#             If True, print detailed debugging information during simulation.
#         rng : optional
#             Random number generator or seed, passed to ``utils._coerce_rng``.
    
#         Returns
#         -------
#         dict
#             Dictionary with the following entries:
    
#             - ``"SteadyStates"`` : list of int  
#               Decimal representations of steady states reached.
    
#             - ``"NumberOfSteadyStates"`` : int  
#               Number of unique steady states found.
    
#             - ``"BasinSizes"`` : list of int  
#               Number of simulations converging to each steady state.
    
#             - ``"TransientTimes"`` : list of list of int  
#               For each steady state, a list of transient lengths (number of update
#               steps before convergence).
    
#             - ``"STGAsynchronous"`` : dict  
#               Partial cache of asynchronous transitions encountered during
#               simulation. Keys are ``(state, node_index)`` and values are successor
#               states (all in decimal form).
    
#             - ``"UpdateQueues"`` : list of list of int  
#               For each simulation, the sequence of visited states (in decimal form).
    
#         Notes
#         -----
#         - Only steady states (fixed points) are detected. If the asynchronous
#           dynamics contain limit cycles, simulations may fail to converge within
#           ``search_depth``.
#         - The returned asynchronous transition graph is incomplete and represents
#           only transitions encountered during sampling.
#         """
#         rng = utils._coerce_rng(rng)
    
#         # --- Initialize initial condition ---
#         if isinstance(initial_condition, int):
#             x0 = utils.dec2bin(initial_condition, self.N)
#             x0dec = initial_condition
#         else:
#             x0 = np.asarray(initial_condition, dtype=int)
#             if x0.shape[0] != self.N:
#                 raise ValueError(
#                     f"Initial condition must have length {self.N}, got {x0.shape[0]}."
#                 )
#             x0dec = utils.bin2dec(x0)
    
#         # --- Handle stochastic weights ---
#         if stochastic_weights is not None:
#             stochastic_weights = np.asarray(stochastic_weights, dtype=float)
#             if stochastic_weights.shape[0] != self.N:
#                 raise ValueError("stochastic_weights must have length N.")
#             if np.any(stochastic_weights <= 0):
#                 raise ValueError("stochastic_weights must be strictly positive.")
#             stochastic_weights = stochastic_weights / stochastic_weights.sum()
    
#         # --- Bookkeeping ---
#         STG_async: dict[tuple[int, int], int] = {}
#         steady_states: list[int] = []
#         basin_sizes: list[int] = []
#         transient_times: list[list[int]] = []
#         steady_state_dict: dict[int, int] = {}
#         queues: list[list[int]] = []
    
#         # --- Simulations ---
#         for iteration in range(n_simulations):
#             x = x0.copy()
#             xdec = x0dec
#             queue = [xdec]
    
#             for step in range(search_depth):
#                 found_new_state = False
    
#                 # If already known steady state, stop
#                 if xdec in steady_state_dict:
#                     idx = steady_state_dict[xdec]
#                     basin_sizes[idx] += 1
#                     transient_times[idx].append(step)
#                     queues.append(queue)
#                     break
    
#                 # Choose update order
#                 if stochastic_weights is None:
#                     update_order = rng.permutation(self.N)
#                 else:
#                     update_order = rng.choice(
#                         self.N, size=self.N, replace=False, p=stochastic_weights
#                     )
    
#                 for i in map(int, update_order):
#                     try:
#                         fxdec = STG_async[(xdec, i)]
#                     except KeyError:
#                         fx_i = self.update_single_node(i, x[self.I[i]])
#                         if fx_i > x[i]:
#                             fxdec = xdec + 2 ** (self.N - 1 - i)
#                             x[i] = 1
#                         elif fx_i < x[i]:
#                             fxdec = xdec - 2 ** (self.N - 1 - i)
#                             x[i] = 0
#                         else:
#                             fxdec = xdec
#                         STG_async[(xdec, i)] = fxdec
    
#                     if fxdec != xdec:
#                         xdec = fxdec
#                         queue.append(xdec)
#                         found_new_state = True
#                         break
    
#                 if debug:
#                     print(iteration, step, i, found_new_state, xdec, x)
    
#                 if not found_new_state:
#                     # New steady state reached
#                     if xdec in steady_state_dict:
#                         idx = steady_state_dict[xdec]
#                         basin_sizes[idx] += 1
#                         transient_times[idx].append(step)
#                     else:
#                         steady_state_dict[xdec] = len(steady_states)
#                         steady_states.append(xdec)
#                         basin_sizes.append(1)
#                         transient_times.append([step])
#                     queues.append(queue)
#                     break
    
#             if debug:
#                 print()
    
#         if sum(basin_sizes) < n_simulations:
#             print(
#                 f"Warning: only {sum(basin_sizes)} of the {n_simulations} simulations "
#                 "reached a steady state. Consider increasing search_depth. "
#                 "The network may contain asynchronous limit cycles."
#             )
    
#         return {
#             "SteadyStates": steady_states,
#             "NumberOfSteadyStates": len(steady_states),
#             "BasinSizes": basin_sizes,
#             "TransientTimes": transient_times,
#             "STGAsynchronous": STG_async,
#             "UpdateQueues": queues,
#         }


    def get_attractors_synchronous(
        self,
        n_simulations: int = 500,
        initial_sample_points: Sequence[int | Sequence[int]] | None = None,
        n_steps_timeout: int = 1000,
        initial_sample_points_are_vectors: bool = False,
        *,
        rng = None,
    ) -> dict:
        """
        Approximate synchronous attractors of a multistate network via sampling.
    
        This method estimates the synchronous attractors (fixed points and cycles)
        of a multistate network by simulating synchronous updates from a collection
        of initial states. For each simulation, the network is updated until an
        attractor is reached or a maximum number of update steps is exceeded.
    
        The method is sampling-based and does *not* guarantee that all attractors
        are found. Basin sizes are lower-bound estimates based on the sampled
        initial conditions.
    
        Parameters
        ----------
        n_simulations : int, optional
            Number of random initial conditions to sample (default is 500). 
            Ignored if ``initial_sample_points`` is provided.
        initial_sample_points : sequence of int or sequence of sequence of int, optional
            Initial states to use. If provided, its length determines the number
            of simulations. Interpretation depends on
            ``initial_sample_points_are_vectors``.
        n_steps_timeout : int, optional
            Maximum number of synchronous update steps per simulation before
            declaring a timeout (default is 1000).
        initial_sample_points_are_vectors : bool, optional
            If True, ``initial_sample_points`` are interpreted as mixed-radix
            vectors; otherwise (default) they are interpreted as decimal-encoded
            states.
        rng : optional
            Random number generator or seed, passed to ``utils._coerce_rng``.
    
        Returns
        -------
        dict
            Dictionary with the following entries:
    
            - ``"Attractors"`` : list of list of int  
              Attractors found, each represented as a list of decimal states
              (cycles are given in cyclic order).
    
            - ``"NumberOfAttractors"`` : int  
              Number of distinct attractors found.
    
            - ``"BasinSizes"`` : list of int  
              Number of sampled initial conditions converging to each attractor.
    
            - ``"AttractorID"`` : dict  
              Mapping from visited states (decimal) to attractor index.
    
            - ``"InitialSamplePoints"`` : list of int  
              Decimal initial states used for sampling.
    
            - ``"STG"`` : dict  
              Sampled synchronous state transition graph
              (state → successor state).
    
            - ``"NumberOfTimeouts"`` : int  
              Number of simulations that did not converge within
              ``n_steps_timeout``.
    
        Notes
        -----
        - This method is intended for networks with long transient dynamics, where
          exhaustive synchronous analysis is infeasible.
        - Basin sizes are *sampling-based estimates* and should not be interpreted
          as exact proportions of the state space.
        - There is no guarantee that all attractors are found. 
        """
        rng = utils._coerce_rng(rng)
    
        # --- Bookkeeping ---
        dictF: dict[int, int] = {}        # memorized synchronous transitions
        attractors: list[list[int]] = []  # attractor cycles
        basin_sizes: list[int] = []       # basin counts
        attr_dict: dict[int, int] = {}    # state -> attractor index
        STG: dict[int, int] = {}           # sampled synchronous STG
        n_timeout = 0
        sampled_points: list[int] = []
    
        initial_sample_points_empty = initial_sample_points is None
        if not initial_sample_points_empty:
            n_simulations = len(initial_sample_points)
    
        # --- Main simulation loop ---
        for sim_idx in range(n_simulations):
            # Initialize state
            if initial_sample_points_empty:
                x = rng.randint(self.R, size = self.N)
                xdec = utils.mix2dec(x, self.R)
                sampled_points.append(xdec)
            else:
                if initial_sample_points_are_vectors:
                    x = np.asarray(initial_sample_points[sim_idx], dtype=np.uint8)
                    if x.shape[0] != self.N:
                        raise ValueError(
                            f"Initial state must have length {self.N}, got {x.shape[0]}."
                        )
                    xdec = utils.mix2dec(x, self.R)
                else:
                    xdec = int(initial_sample_points[sim_idx])
                    x = np.array(utils.dec2bin(xdec, self.N), dtype=np.uint8)
    
            visited = {xdec: 0}
            trajectory = [xdec]
            count = 0
    
            # --- Iterate until attractor or timeout ---
            while count < n_steps_timeout:
                if xdec in dictF:
                    fxdec = dictF[xdec]
                else:
                    fx = self._update_network_synchronously_unchecked(x)
    
                    fxdec = utils.bin2dec(fx)
                    dictF[xdec] = fxdec
                    x = fx
    
                # record sampled STG edge (first visit only)
                if count == 0:
                    STG[xdec] = fxdec
    
                # already assigned to known attractor
                if fxdec in attr_dict:
                    idx_attr = attr_dict[fxdec]
                    basin_sizes[idx_attr] += 1
                    for s in trajectory:
                        attr_dict[s] = idx_attr
                    break
    
                # new attractor detected
                if fxdec in visited:
                    cycle_start = visited[fxdec]
                    attractor_states = trajectory[cycle_start:]
                    attractors.append(attractor_states)
                    basin_sizes.append(1)
                    idx_attr = len(attractors) - 1
                    for s in attractor_states:
                        attr_dict[s] = idx_attr
                    break
    
                # continue traversal
                visited[fxdec] = len(trajectory)
                trajectory.append(fxdec)
                xdec = fxdec
                count += 1
    
                if count == n_steps_timeout:
                    n_timeout += 1
                    break
    
        return {
            "Attractors": attractors,
            "NumberOfAttractors": len(attractors),
            "BasinSizes": basin_sizes,
            "AttractorID": attr_dict,
            "InitialSamplePoints": (
                sampled_points if initial_sample_points_empty else list(initial_sample_points)
            ),
            "STG": STG,
            "NumberOfTimeouts": n_timeout,
        }


    def compute_synchronous_state_transition_graph(self) -> None:
        """
        Compute the exact synchronous state transition graph (STG)
        """
        # this is a slow implementation of this functionality
        # to be improved in the future
        states_dec = list(range(self.R.prod(0, int)))
        next_states = np.zeros_like(states_dec, np.uint8)
        for i in states_dec:
            next_states[i] = self.update_network_synchronously(utils.dec2mix(i, self.R))
        self.STG = next_states
#    
#         states = utils.get_left_side_of_truth_table(self.N)
    
#         # 2. Allocate next-state matrix
#         next_states = np.zeros_like(states, dtype=np.uint8)
    
#         # Binary-to-decimal weights
#         powers_of_two = (1 << np.arange(self.N))[::-1]
    
#         # 3. Compute next state for each node
#         for j, bf in enumerate(self.F):
#             regulators = self.I[j]
    
#             if len(regulators) == 0:
#                 # Constant node
#                 next_states[:, j] = bf.f[0]
#                 continue
    
#             subspace = states[:, regulators]
#             idx = np.dot(subspace, powers_of_two[-len(regulators):])
#             next_states[:, j] = bf.f[idx]
    
#         # 4. Convert next states to decimal
#         self.STG = np.dot(next_states, powers_of_two).astype(np.int64)


    def get_attractors_synchronous_exact(
        self,
    ) -> dict:
        """
        Compute all attractors and their exact basin sizes under synchronous updating.
    
        This method computes the exact synchronous state transition graph (STG) and
        analyzes it as a functional graph on all states. All attractors (cycles),
        their basin sizes, and the attractor reached from each state are determined
        exactly.
    
        This computation requires memory and time proportional to the network size and is
        intended for small-to-moderate networks only.
    
        Returns
        -------
        dict
            Dictionary with keys:
    
            - Attractors : list[list[int]]
                Each attractor represented as a list of decimal states forming a cycle.
            - NumberOfAttractors : int
                Total number of attractors.
            - BasinSizes : np.ndarray[float]
                Fraction of all states belonging to each attractor basin.
            - AttractorID : np.ndarray[int]
                For each of the ``2**N`` states, the index of the attractor it reaches.
            - STG : np.ndarray[int]
                The synchronous state transition graph.
        """
        if self.STG is None:
            self.compute_synchronous_state_transition_graph()
    
        attractors = []
        attractor_id = -np.ones(self.R.prod(0, int), dtype=np.int32)
        basin_sizes = []
        n_attr = 0

        for xdec in range(self.R.prod(0, int)):
            if attractor_id[xdec] != -1:
                continue

            cur = xdec
            queue = [cur]

            while True:
                fxdec = int(self.STG[cur])

                if attractor_id[fxdec] != -1:
                    idx_attr = attractor_id[fxdec]
                    basin_sizes[idx_attr] += len(queue)
                    for q in queue:
                        attractor_id[q] = idx_attr
                    break

                if fxdec in queue:
                    idx = queue.index(fxdec)
                    cycle = queue[idx:]
                    attractors.append(cycle)
                    basin_sizes.append(len(queue))
                    for q in queue:
                        attractor_id[q] = n_attr
                    n_attr += 1
                    break

                queue.append(fxdec)
                cur = fxdec
    
        basin_sizes = np.array(basin_sizes, dtype=np.float64) / self.R.prod(0, int)
    
        return {
            "Attractors": attractors,
            "NumberOfAttractors": len(attractors),
            "BasinSizes": basin_sizes,
            "AttractorID": attractor_id,
            "STG": self.STG,
        }

    
#     def get_transient_lengths_exact(
#         self,
#         use_numba : bool = True
#     ) -> np.ndarray:
#         """
#         Compute exact transient length using:
#           - Full STG from get_attractors_synchronous_exact()
#           - Attractors (cycle states) from get_attractors_synchronous_exact()
    
#         This avoids indegree-pruning because cycle states are given explicitly.
#         """
#         attractor_info = self.get_attractors_synchronous_exact(use_numba=use_numba)

#         stg = self.STG                              # full mapping: successor(s)
#         attractors = attractor_info["Attractors"]   # list of cycles
        
#         if __LOADED_NUMBA__ and use_numba:
#             is_attr_mask = np.full(2**self.N, 0, dtype=np.uint8)
        
#             for i, states in enumerate(attractors):
#                 states_arr = np.asarray(states, dtype=np.int64)
#                 is_attr_mask[states_arr] = 1
#             return _transient_lengths_functional_numba(
#                 self.STG.astype(np.int64, copy=False),
#                 is_attr_mask
#             )
        
        
#         # Normalize STG to an integer array/list succ where succ[u] = v
#         if isinstance(stg, np.ndarray):
#             succ = stg.astype(int, copy=False)
#             n = int(succ.shape[0])
#         else:
#             succ = list(stg)
#             n = len(succ)
    
#         # Build reverse adjacency list rev[v] = all u such that u -> v
#         rev = [[] for _ in range(n)]
#         for u in range(n):
#             v = int(succ[u])
#             if v < 0 or v >= n:
#                 raise ValueError(f"Invalid successor: {u} -> {v}")
#             rev[v].append(u)
    
#         # Initialize distances: all cycle states have transient length 0
#         dist = np.full(n, -1, dtype=np.int64)
#         bfs = deque()
    
#         for cycle in attractors:
#             for s in cycle:
#                 if dist[s] == -1:
#                     dist[s] = 0
#                     bfs.append(s)
    
#         # Multi-source BFS outward from cycle states
#         while bfs:
#             v = bfs.popleft()
#             for u in rev[v]:
#                 if dist[u] == -1:
#                     dist[u] = dist[v] + 1
#                     bfs.append(u)
    
#         # If STG is complete, every state must get a distance
#         if any(d < 0 for d in dist):
#             raise RuntimeError("Some states did not receive a transient length. Is STG complete?")
        
#         return np.array(dist,dtype=int)

#     ## Robustness measures: synchronous Derrida value, entropy of basin size distribution, coherence, fragility
#     def get_attractors_and_robustness_synchronous_exact(
#         self, 
#         use_numba: bool = True,
#         get_stratified_coherences : bool = False
#     ) -> dict:
#         """
#         Compute attractors and exact robustness measures of a synchronously
#         updated Boolean network.

#         This method constructs the exact synchronous state transition graph
#         (STG) on ``2**N`` states and analyzes it as a functional graph. All
#         attractors (cycles), basin sizes, and the attractor reached from each
#         state are determined exactly. Based on this decomposition, exact
#         coherence and fragility measures are computed for the full network,
#         for each basin of attraction, and for each attractor.

#         Optionally, coherence can be stratified by the transient length
#         (distance from the attractor) of each state, allowing robustness to be
#         analyzed as a function of how far states lie from their eventual
#         attractor.

#         This computation requires memory and time proportional to ``2**N`` and
#         is intended for small-to-moderate networks. When Numba is enabled,
#         exact and stratified robustness measures remain feasible up to
#         moderate values of ``N`` (e.g., ``N ≈ 20`` on typical hardware).

#         Parameters
#         ----------
#         use_numba : bool, optional
#             If True (default) and Numba is available, compiled kernels are used
#             for robustness and transient-length computations, resulting in
#             substantial speedups.
#         get_stratified_coherences : bool, optional
#             If True, coherence is additionally computed as a function of the
#             transient length (distance to the attractor) of each state.
#             When Numba is enabled, this option incurs only modest additional
#             computational cost. Default is False.

#         Returns
#         -------
#         dict
#             Dictionary with the following keys:

#             - Attractors : list[list[int]]
#                 Each attractor represented as a list of decimal states forming
#                 a cycle.
#             - NumberOfAttractors : int
#                 Total number of attractors.
#             - BasinSizes : np.ndarray of float
#                 Fraction of all states belonging to each attractor basin.
#             - AttractorID : np.ndarray of int
#                 For each of the ``2**N`` states, the index of the attractor it
#                 eventually reaches.
#             - Coherence : float
#                 Exact global network coherence.
#             - Fragility : float
#                 Exact global network fragility.
#             - BasinCoherence : np.ndarray of float
#                 Exact coherence of each basin of attraction.
#             - BasinFragility : np.ndarray of float
#                 Exact fragility of each basin of attraction.
#             - AttractorCoherence : np.ndarray of float
#                 Exact coherence of each attractor.
#             - AttractorFragility : np.ndarray of float
#                 Exact fragility of each attractor.

#             If ``get_stratified_coherences`` is True, the dictionary additionally
#             contains:

#             - StratifiedCoherences : np.ndarray of float
#                 Coherence values stratified by attractor and transient length.
#             - DistanceFromAttractorCount : np.ndarray of int
#                 Number of state–hypercube-edge incidences contributing to each
#                 stratified coherence entry.
#             - DistanceFromAttractor : np.ndarray of int
#                 Transient length (distance to attractor) for each state.
#         """
    
#         # ------------------------------------------------------------------
#         # 0) Attractors and basins
#         # ------------------------------------------------------------------
#         result = self.get_attractors_synchronous_exact(use_numba=use_numba)
    
#         attractors = result["Attractors"]
#         n_attractors = int(result["NumberOfAttractors"])
    
#         basin_sizes = np.asarray(result["BasinSizes"], dtype=np.float64)
#         attractor_id = np.asarray(result["AttractorID"], dtype=np.int64)
    
#         n_states = 1 << self.N
    
#         # ------------------------------------------------------------------
#         # Single-attractor shortcut
#         # ------------------------------------------------------------------
#         if n_attractors == 1:
#             return {
#                 "Attractors": attractors,
#                 "NumberOfAttractors": 1,
#                 "BasinSizes": basin_sizes,
#                 "AttractorID": attractor_id,
#                 "Coherence": 1.0,
#                 "Fragility": 0.0,
#                 "BasinCoherence": np.ones(1, dtype=np.float64),
#                 "BasinFragility": np.zeros(1, dtype=np.float64),
#                 "AttractorCoherence": np.ones(1, dtype=np.float64),
#                 "AttractorFragility": np.zeros(1, dtype=np.float64),
#             }
    
#         # ------------------------------------------------------------------
#         # 1) Attractor membership and lengths
#         # ------------------------------------------------------------------
#         is_attr_mask = np.zeros(n_states, dtype=np.uint8)
#         len_attractors = np.empty(n_attractors, dtype=np.int64)
    
#         for i, states in enumerate(attractors):
#             states_arr = np.asarray(states, dtype=np.int64)
#             len_attractors[i] = states_arr.size
#             is_attr_mask[states_arr] = 1
    
#         # ------------------------------------------------------------------
#         # 2) Mean binary vector per attractor
#         # ------------------------------------------------------------------
#         mean_states_attractors = np.empty((n_attractors, self.N), dtype=np.float64)
    
#         for i, states in enumerate(attractors):
#             if len(states) == 1:
#                 mean_states_attractors[i] = np.asarray(
#                     utils.dec2bin(states[0], self.N), dtype=np.float64
#                 )
#             else:
#                 arr = np.asarray(
#                     [utils.dec2bin(s, self.N) for s in states], dtype=np.float64
#                 )
#                 mean_states_attractors[i] = arr.mean(axis=0)
    
#         # ------------------------------------------------------------------
#         # 3) Distance matrix between attractors
#         # ------------------------------------------------------------------
#         diff = mean_states_attractors[:, None, :] - mean_states_attractors[None, :, :]
#         distance_between_attractors = np.sum(np.abs(diff), axis=2)
#         distance_between_attractors = np.asarray(
#             distance_between_attractors / float(self.N), dtype=np.float64
#         )
    
#         # ------------------------------------------------------------------
#         # 4) Hypercube edge traversal
#         # ------------------------------------------------------------------
#         if __LOADED_NUMBA__ and use_numba:
#             if get_stratified_coherences:
#                 distances_from_attractor = _transient_lengths_functional_numba(
#                     self.STG.astype(np.int64, copy=False),
#                     is_attr_mask
#                 )
#                 max_distance_from_attractor = int(distances_from_attractor.max())
        
#                 (
#                     basin_coherences,
#                     basin_fragilities,
#                     attractor_coherences,
#                     attractor_fragilities,
#                     stratified_coherences,
#                     n_states_with_specific_distance_from_attractor,
#                 ) = _robustness_edge_traversal_numba_stratified(
#                     int(self.N),
#                     attractor_id,
#                     is_attr_mask,
#                     distance_between_attractors,
#                     distances_from_attractor,
#                     max_distance_from_attractor,
#                 )
                    
#                 stratified_coherences = np.asarray(stratified_coherences, dtype=np.float64)
#                 n_states_with_specific_distance_from_attractor = np.asarray(n_states_with_specific_distance_from_attractor, dtype=int)
#             else:
#                 (
#                     basin_coherences,
#                     basin_fragilities,
#                     attractor_coherences,
#                     attractor_fragilities,
#                 ) = _robustness_edge_traversal_numba(
#                     int(self.N),
#                     attractor_id,
#                     is_attr_mask,
#                     distance_between_attractors,
#                 )
    
#             basin_coherences = np.asarray(basin_coherences, dtype=np.float64)
#             basin_fragilities = np.asarray(basin_fragilities, dtype=np.float64)
#             attractor_coherences = np.asarray(attractor_coherences, dtype=np.float64)
#             attractor_fragilities = np.asarray(attractor_fragilities, dtype=np.float64)
    
#         else:
#             basin_coherences = np.zeros(n_attractors, dtype=np.float64)
#             basin_fragilities = np.zeros(n_attractors, dtype=np.float64)
#             attractor_coherences = np.zeros(n_attractors, dtype=np.float64)
#             attractor_fragilities = np.zeros(n_attractors, dtype=np.float64)
            
#             if get_stratified_coherences:
#                 distances_from_attractor = self.get_transient_lengths_exact(result)
#                 max_distance_from_attractor = max(distances_from_attractor)
#                 stratified_coherences = np.zeros((n_attractors,max_distance_from_attractor+1), dtype=np.float64)
#                 n_states_with_specific_distance_from_attractor = np.zeros((n_attractors,max_distance_from_attractor+1), dtype=int)
                
#             for xdec in range(n_states):
#                 for bitpos in range(self.N):
#                     if (xdec >> bitpos) & 1:
#                         continue
    
#                     ydec = xdec | (1 << bitpos)
    
#                     idx_x = attractor_id[xdec]
#                     idx_y = attractor_id[ydec]
                    
#                     if get_stratified_coherences:
#                         n_states_with_specific_distance_from_attractor[idx_x,distances_from_attractor[xdec]] += 1
#                         n_states_with_specific_distance_from_attractor[idx_y,distances_from_attractor[ydec]] += 1
                        
#                     if idx_x == idx_y:
#                         basin_coherences[idx_x] += 2.0
#                         if is_attr_mask[xdec]:
#                             attractor_coherences[idx_x] += 1.0
#                         if is_attr_mask[ydec]:
#                             attractor_coherences[idx_y] += 1.0
#                         if get_stratified_coherences:
#                             stratified_coherences[idx_x,distances_from_attractor[xdec]] += 1.0
#                             stratified_coherences[idx_y,distances_from_attractor[ydec]] += 1.0
#                     else:
#                         dxy = float(distance_between_attractors[idx_x, idx_y])
#                         basin_fragilities[idx_x] += dxy
#                         basin_fragilities[idx_y] += dxy
#                         if is_attr_mask[xdec]:
#                             attractor_fragilities[idx_x] += dxy
#                         if is_attr_mask[ydec]:
#                             attractor_fragilities[idx_y] += dxy
    
#         # ------------------------------------------------------------------
#         # 5) Normalization
#         # ------------------------------------------------------------------
#         basin_counts = basin_sizes * float(n_states)
    
#         if get_stratified_coherences:
#             n_states_with_specific_distance_from_attractor //= self.N
    
#         for i in range(n_attractors):
#             if basin_counts[i] > 0.0:
#                 basin_coherences[i] /= basin_counts[i] * self.N
#                 basin_fragilities[i] /= basin_counts[i] * self.N
    
#             if len_attractors[i] > 0:
#                 attractor_coherences[i] /= len_attractors[i] * self.N
#                 attractor_fragilities[i] /= len_attractors[i] * self.N
                
#             if get_stratified_coherences:
#                 for d in range(max_distance_from_attractor+1):
#                     if n_states_with_specific_distance_from_attractor[i,d] > 0.0:
#                         stratified_coherences[i,d] /= n_states_with_specific_distance_from_attractor[i,d] * self.N
#                     else:
#                         stratified_coherences[i,d] = np.nan
                        
#         coherence = float(np.dot(basin_sizes, basin_coherences))
#         fragility = float(np.dot(basin_sizes, basin_fragilities))
    
#         # ------------------------------------------------------------------
#         # Final return
#         # ------------------------------------------------------------------
#         return_dict =  {
#             "Attractors": attractors,
#             "NumberOfAttractors": int(n_attractors),
#             "BasinSizes": basin_sizes,
#             "AttractorID": attractor_id,
#             "Coherence": coherence,
#             "Fragility": fragility,
#             "BasinCoherence": basin_coherences,
#             "BasinFragility": basin_fragilities,
#             "AttractorCoherence": attractor_coherences,
#             "AttractorFragility": attractor_fragilities,
#         }
#         if get_stratified_coherences:
#             return_dict['StratifiedCoherences'] = stratified_coherences
#             return_dict['DistanceFromAttractorCount'] = n_states_with_specific_distance_from_attractor
#             return_dict['DistanceFromAttractor'] = distances_from_attractor
#         return return_dict


#     def get_attractors_and_robustness_synchronous(
#         self,
#         n_simulations: int = 500,
#         return_attractor_coherence: bool = True,
#         *,
#         rng=None,
#     ) -> dict:
#         """
#         Approximate attractors and robustness measures under synchronous updating.
    
#         This method samples the attractor landscape by simulating the network from
#         multiple random initial conditions (ICs) and their single-bit perturbations.
#         It returns Monte-Carlo approximations of global coherence, fragility, and a
#         final Hamming-distance-based measure, along with per-basin approximations.
#         Optionally, it additionally estimates attractor-level coherence and fragility
#         by perturbing attractor states found during sampling.
    
#         Notes
#         -----
#         - The attractor set returned is a *lower bound* on the true number of
#           attractors, because only the sampled portion of state space is explored.
#         - For ``N >= 64``, decimal encoding of states may exceed ``np.int64`` and
#           this method uses bitstrings (type ``str``) as state identifiers.
    
#         Parameters
#         ----------
#         n_simulations : int, optional
#             Number of random initial conditions to sample (default is 500). For each
#             IC, the method also simulates one randomly chosen single-bit perturbation.
#         return_attractor_coherence : bool, optional
#             If True (default), also compute attractor-level coherence and fragility
#             by perturbing attractor states found during sampling.
#         rng : None or numpy.random.Generator, optional
#             Random number generator or seed-like object. Passed to
#             ``utils._coerce_rng``.
    
#         Returns
#         -------
#         dict
#             Dictionary with keys:
    
#             - Attractors : list[list[int]] or list[list[str]]
#                 List of discovered attractors, each represented as a list of states
#                 forming a cycle. States are decimals (``int``) for ``N < 64`` and
#                 bitstrings (``str``) for ``N >= 64``.
#             - LowerBoundOfNumberOfAttractors : int
#                 Number of distinct attractors discovered (a lower bound on the true
#                 number of attractors).
#             - BasinSizesApproximation : np.ndarray[float]
#                 Approximate basin size (fraction of sampled trajectories that end in
#                 each attractor). Sums to ~1 over discovered attractors.
#             - CoherenceApproximation : float
#                 Approximate global coherence: probability that a random IC and its
#                 single-bit perturbation reach the same attractor.
#             - FragilityApproximation : float
#                 Approximate global fragility: expected normalized difference between
#                 reached attractors when the IC and perturbation reach different
#                 attractors. Normalized by ``N``.
#             - FinalHammingDistanceApproximation : float
#                 Approximate final Hamming distance between the two periodic
#                 trajectories when comparing the IC and its perturbation. This is a
#                 *distance* in [0, 1], where 0 means identical and 1 means completely
#                 different.
#             - BasinCoherenceApproximation : np.ndarray[float]
#                 Approximate coherence per basin (same definition as coherence but
#                 conditioned on having reached that basin).
#             - BasinFragilityApproximation : np.ndarray[float]
#                 Approximate fragility per basin (same definition as fragility but
#                 conditioned on having reached that basin).
#             - AttractorCoherence : np.ndarray[float], optional
#                 If ``return_attractor_coherence`` is True: estimated attractor-level
#                 coherence (probability that a single-bit perturbation of an attractor
#                 state returns to the same attractor).
#             - AttractorFragility : np.ndarray[float], optional
#                 If ``return_attractor_coherence`` is True: estimated attractor-level
#                 fragility based on differences between the original attractor and the
#                 attractor reached after perturbation.
    
#         References
#         ----------
#         Park, K. H., Costa, F. X., Rocha, L. M., Albert, R., & Rozum, J. C. (2023).
#         Models of cell processes are far from the edge of chaos. PRX Life, 1(2), 023009.
    
#         Bavisetty, V. S. N., Wheeler, M., & Kadelka, C. (2025).
#         Attractors are less stable than their basins: Canalization creates a coherence
#         gap in gene regulatory networks. bioRxiv 2025-11.
#         """
#         rng = utils._coerce_rng(rng)
    
#         def lcm(a: int, b: int) -> int:
#             return abs(a * b) // math.gcd(a, b)
    
#         # ------------------------------------------------------------------
#         # Initialization
#         # ------------------------------------------------------------------
#         dictF = {}
#         attractors = []
#         ICs_per_attractor_state = []
#         basin_sizes = []
#         attractor_dict = {}
#         attractor_state_dict = []
#         distance_from_attractor_state_dict = []
#         counter_phase_shifts = []
    
#         powers_of_2s = [
#             np.asarray([2**i for i in range(NN)][::-1], dtype=np.int64)
#             for NN in range(max(self.indegrees) + 1)
#         ]
    
#         if self.N < 64:
#             powers_of_2 = np.asarray([2**i for i in range(self.N)][::-1], dtype=np.int64)
    
#         robustness_approximation = 0
#         fragility_sum = 0.0
#         basin_robustness = defaultdict(float)
#         basin_fragility = defaultdict(float)
#         final_hamming_distance_approximation = 0.0
    
#         mean_states_attractors = []
#         states_attractors = []
    
#         # ------------------------------------------------------------------
#         # Sampling phase
#         # ------------------------------------------------------------------
#         for _ in range(n_simulations):
#             index_attractors = []
#             index_within_attr = []
#             dist_from_attr = []
    
#             for j in range(2):
#                 if j == 0:
#                     x = rng.integers(2, size=self.N, dtype=np.uint8)
#                     if self.N < 64:
#                         xdec = int(np.dot(x, powers_of_2))
#                     else:
#                         xdec = "".join(str(int(b)) for b in x)
#                     x_old = x.copy()
#                 else:
#                     x = x_old.copy()
#                     bit = int(rng.integers(self.N))
#                     x[bit] ^= 1
#                     if self.N < 64:
#                         xdec = int(np.dot(x, powers_of_2))
#                     else:
#                         xdec = "".join(str(int(b)) for b in x)
    
#                 queue = [xdec]
    
#                 try:
#                     idx_attr = attractor_dict[xdec]
#                 except KeyError:
#                     while True:
#                         try:
#                             fxdec = dictF[xdec]
#                         except KeyError:
#                             fx = np.empty(self.N, dtype=np.uint8)
#                             for jj in range(self.N):
#                                 if self.indegrees[jj] > 0:
#                                     fx[jj] = self.F[jj].f[
#                                         int(
#                                             np.dot(
#                                                 x[self.I[jj]],
#                                                 powers_of_2s[self.indegrees[jj]],
#                                             )
#                                         )
#                                     ]
#                                 else:
#                                     fx[jj] = self.F[jj].f[0]
    
#                             if self.N < 64:
#                                 fxdec = int(np.dot(fx, powers_of_2))
#                             else:
#                                 fxdec = "".join(str(int(b)) for b in fx)
    
#                             dictF[xdec] = fxdec
    
#                         try:
#                             idx_attr = attractor_dict[fxdec]
#                             idx_state = attractor_state_dict[idx_attr][fxdec]
#                             dist_state = distance_from_attractor_state_dict[idx_attr][fxdec]
    
#                             attractor_dict.update({q: idx_attr for q in queue})
#                             attractor_state_dict[idx_attr].update(
#                                 {q: idx_state for q in queue}
#                             )
#                             distance_from_attractor_state_dict[idx_attr].update(
#                                 {
#                                     q: d
#                                     for q, d in zip(
#                                         queue,
#                                         range(len(queue) + dist_state, dist_state, -1),
#                                     )
#                                 }
#                             )
#                             break
    
#                         except KeyError:
#                             if fxdec in queue:
#                                 idx = queue.index(fxdec)
#                                 idx_attr = len(attractors)
    
#                                 attractors.append(queue[idx:])
#                                 basin_sizes.append(1)
#                                 ICs_per_attractor_state.append(
#                                     [0] * len(attractors[-1])
#                                 )
#                                 counter_phase_shifts.append(
#                                     [0] * len(attractors[-1])
#                                 )
    
#                                 attractor_dict.update({q: idx_attr for q in queue})
#                                 attractor_state_dict.append(
#                                     {
#                                         q: (0 if q in queue[:idx] else queue[idx:].index(q))
#                                         for q in queue
#                                     }
#                                 )
#                                 distance_from_attractor_state_dict.append(
#                                     {
#                                         q: (idx - queue.index(q))
#                                         if q in queue[:idx]
#                                         else 0
#                                         for q in queue
#                                     }
#                                 )
    
#                                 if len(attractors[-1]) == 1:
#                                     fp = (
#                                         np.asarray(
#                                             utils.dec2bin(queue[idx], self.N),
#                                             dtype=np.float64,
#                                         )
#                                         if self.N < 64
#                                         else np.asarray(list(queue[idx]), dtype=np.float64)
#                                     )
#                                     states_attractors.append(fp.reshape(1, self.N))
#                                     mean_states_attractors.append(fp)
#                                 else:
#                                     lc = (
#                                         np.asarray(
#                                             [
#                                                 utils.dec2bin(s, self.N)
#                                                 for s in queue[idx:]
#                                             ],
#                                             dtype=np.float64,
#                                         )
#                                         if self.N < 64
#                                         else np.asarray(
#                                             [list(s) for s in queue[idx:]],
#                                             dtype=np.float64,
#                                         )
#                                     )
#                                     states_attractors.append(lc)
#                                     mean_states_attractors.append(lc.mean(axis=0))
#                                 break
#                             else:
#                                 x = fx.copy()
#                                 queue.append(fxdec)
#                                 xdec = fxdec
    
#                 index_attractors.append(idx_attr)
#                 index_within_attr.append(attractor_state_dict[idx_attr][xdec])
#                 dist_from_attr.append(
#                     distance_from_attractor_state_dict[idx_attr][xdec]
#                 )
    
#                 basin_sizes[idx_attr] += 1
#                 ICs_per_attractor_state[idx_attr][
#                     attractor_state_dict[idx_attr][xdec]
#                 ] += 1
    
#             if index_attractors[0] == index_attractors[1]:
#                 robustness_approximation += 1
#                 basin_robustness[index_attractors[0]] += 1
#                 ps = max(index_within_attr) - min(index_within_attr)
#                 counter_phase_shifts[index_attractors[0]][ps] += 1
#             else:
#                 d = np.sum(
#                     np.abs(
#                         mean_states_attractors[index_attractors[0]]
#                         - mean_states_attractors[index_attractors[1]]
#                     )
#                 )
#                 fragility_sum += d
#                 basin_fragility[index_attractors[0]] += d
    
#                 L = lcm(
#                     len(attractors[index_attractors[0]]),
#                     len(attractors[index_attractors[1]]),
#                 )
    
#                 s0 = states_attractors[index_attractors[0]]
#                 s1 = states_attractors[index_attractors[1]]
    
#                 p0 = np.tile(s0, (L // len(s0) + 1, 1))[
#                     index_within_attr[0] : index_within_attr[0] + L
#                 ]
#                 p1 = np.tile(s1, (L // len(s1) + 1, 1))[
#                     index_within_attr[1] : index_within_attr[1] + L
#                 ]
    
#                 final_hamming_distance_approximation += np.mean(p0 == p1)
    
#         # ------------------------------------------------------------------
#         # Aggregation
#         # ------------------------------------------------------------------
#         lower_bound_number_of_attractors = len(attractors)
    
#         approximate_basin_sizes = (
#             np.asarray(basin_sizes, dtype=np.float64)
#             / (2.0 * float(n_simulations))
#         )
    
#         approximate_coherence = robustness_approximation / float(n_simulations)
#         approximate_fragility = fragility_sum / float(n_simulations) / float(self.N)
    
#         approximate_basin_coherence = np.asarray(
#             [
#                 2.0 * basin_robustness[i] / basin_sizes[i]
#                 for i in range(lower_bound_number_of_attractors)
#             ],
#             dtype=np.float64,
#         )
    
#         approximate_basin_fragility = np.asarray(
#             [
#                 2.0 * basin_fragility[i] / basin_sizes[i] / float(self.N)
#                 for i in range(lower_bound_number_of_attractors)
#             ],
#             dtype=np.float64,
#         )
    
#         final_hamming_distance_approximation /= float(n_simulations)
    
#         results = [
#             attractors,
#             lower_bound_number_of_attractors,
#             approximate_basin_sizes,
#             approximate_coherence,
#             approximate_fragility,
#             final_hamming_distance_approximation,
#             approximate_basin_coherence,
#             approximate_basin_fragility,
#         ]
    
#         if not return_attractor_coherence:
#             return dict(
#                 zip(
#                     [
#                         "Attractors",
#                         "LowerBoundOfNumberOfAttractors",
#                         "BasinSizesApproximation",
#                         "CoherenceApproximation",
#                         "FragilityApproximation",
#                         "FinalHammingDistanceApproximation",
#                         "BasinCoherenceApproximation",
#                         "BasinFragilityApproximation",
#                     ],
#                     results,
#                 )
#             )
    
#         # ------------------------------------------------------------------
#         # Attractor-level coherence / fragility
#         # ------------------------------------------------------------------
#         attractor_coherence = np.zeros(lower_bound_number_of_attractors, dtype=np.float64)
#         attractor_fragility = np.zeros(lower_bound_number_of_attractors, dtype=np.float64)
    
#         attractors_original = attractors[:]
    
#         for idx0, attractor in enumerate(attractors_original):
#             for state in attractor:
#                 for i in range(self.N):
#                     x = (
#                         np.asarray(utils.dec2bin(state, self.N), dtype=np.uint8)
#                         if self.N < 64
#                         else np.asarray(list(state), dtype=np.uint8)
#                     )
#                     x[i] ^= 1
    
#                     if self.N < 64:
#                         xdec = int(np.dot(x, powers_of_2))
#                     else:
#                         xdec = "".join(str(int(b)) for b in x)
    
#                     try:
#                         idx1 = attractor_dict[xdec]
#                     except KeyError:
#                         # --- safe forward-walk without touching basin counts
#                         queue = [xdec]
#                         x_local = x.copy()
#                         while True:
#                             try:
#                                 fxdec = dictF[xdec]
#                             except KeyError:
#                                 fx = np.empty(self.N, dtype=np.uint8)
#                                 for jj in range(self.N):
#                                     if self.indegrees[jj] > 0:
#                                         fx[jj] = self.F[jj].f[
#                                             int(
#                                                 np.dot(
#                                                     x_local[self.I[jj]],
#                                                     powers_of_2s[self.indegrees[jj]],
#                                                 )
#                                             )
#                                         ]
#                                     else:
#                                         fx[jj] = self.F[jj].f[0]
    
#                                 if self.N < 64:
#                                     fxdec = int(np.dot(fx, powers_of_2))
#                                 else:
#                                     fxdec = "".join(str(int(b)) for b in fx)
    
#                                 dictF[xdec] = fxdec
    
#                             if fxdec in attractor_dict:
#                                 idx1 = attractor_dict[fxdec]
#                                 break
    
#                             if fxdec in queue:
#                                 idx1 = len(attractors)
#                                 attractors.append(queue[queue.index(fxdec):])
#                                 attractor_dict.update(
#                                     {q: idx1 for q in queue}
#                                 )
#                                 break
    
#                             queue.append(fxdec)
#                             xdec = fxdec
#                             x_local = fx.copy()
    
#                     if idx0 == idx1:
#                         attractor_coherence[idx0] += 1.0
#                     else:
#                         attractor_fragility[idx0] += np.sum(
#                             np.abs(
#                                 mean_states_attractors[idx0]
#                                 - mean_states_attractors[idx1]
#                             )
#                         )
    
#         attractor_coherence /= (
#             float(self.N)
#             * np.asarray(list(map(len, attractors_original)), dtype=np.float64)
#         )
    
#         attractor_fragility /= (
#             float(self.N) ** 2
#             * np.asarray(list(map(len, attractors_original)), dtype=np.float64)
#         )
    
#         results[0] = attractors_original
    
#         return dict(
#             zip(
#                 [
#                     "Attractors",
#                     "LowerBoundOfNumberOfAttractors",
#                     "BasinSizesApproximation",
#                     "CoherenceApproximation",
#                     "FragilityApproximation",
#                     "FinalHammingDistanceApproximation",
#                     "BasinCoherenceApproximation",
#                     "BasinFragilityApproximation",
#                     "AttractorCoherence",
#                     "AttractorFragility",
#                 ],
#                 results + [attractor_coherence, attractor_fragility],
#             )
#         )
    
    
#     def get_derrida_value(
#         self,
#         n_simulations: int = 1000,
#         exact: bool = False,
#         use_numba: bool = True,
#         *,
#         rng=None,
#     ) -> float:
#         """
#         Compute the Derrida value of a Boolean network.
    
#         The Derrida value measures the average Hamming distance between the
#         one-step synchronous updates of two states that differ by a single-bit
#         perturbation. It quantifies the short-term sensitivity of the network
#         dynamics to small perturbations.
    
#         If ``exact`` is True, the Derrida value is computed exactly as the mean
#         (unnormalized) average sensitivity of the Boolean update functions.
#         Otherwise, it is approximated via Monte Carlo simulation.
    
#         Parameters
#         ----------
#         n_simulations : int, optional
#             Number of Monte Carlo simulations to perform (default is 1000).
#             Ignored if ``exact`` is True.
#         exact : bool, optional
#             If True, compute the exact Derrida value. If False (default),
#             approximate the Derrida value using Monte Carlo simulation.
#         use_numba : bool, optional
#             If True (default) and Numba is available, use a compiled kernel for
#             Monte Carlo simulation.
#         rng : None or np.random.Generator, optional
#             Random number generator, passed through ``utils._coerce_rng``.
    
#         Returns
#         -------
#         float
#             The Derrida value, defined as the average Hamming distance after
#             one synchronous update following a single-bit perturbation.
    
#         References
#         ----------
#         Derrida, B., & Pomeau, Y. (1986).
#         Random networks of automata: a simple annealed approximation.
#         *Europhysics Letters*, 1(2), 45.
#         """
    
#         # ------------------------------------------------------------------
#         # Exact computation
#         # ------------------------------------------------------------------
#         if exact:
#             return float(
#                 np.mean(
#                     [
#                         bf.get_average_sensitivity(
#                             exact=True, normalized=False
#                         )
#                         for bf in self.F
#                     ]
#                 )
#             )
    
#         # ------------------------------------------------------------------
#         # Monte Carlo approximation
#         # ------------------------------------------------------------------
#         rng = utils._coerce_rng(rng)
    
#         if __LOADED_NUMBA__ and use_numba:
#             # Prepare Numba-friendly inputs
#             F_array_list = List(
#                 [np.asarray(bf.f, dtype=np.uint8) for bf in self.F]
#             )
#             I_array_list = List(
#                 [np.asarray(regs, dtype=np.int64) for regs in self.I]
#             )
    
#             seed = int(rng.integers(0, 2**31 - 1))
    
#             return float(
#                 _derrida_simulation(
#                     F_array_list,
#                     I_array_list,
#                     int(self.N),
#                     int(n_simulations),
#                     seed,
#                 )
#             )
    
#         # ------------------------------------------------------------------
#         # Pure Python fallback
#         # ------------------------------------------------------------------
#         total_dist: float = 0.0
    
#         for _ in range(int(n_simulations)):
#             x = rng.integers(0, 2, size=self.N, dtype=np.uint8)
#             y = x.copy()
    
#             idx = int(rng.integers(0, self.N))
#             y[idx] ^= np.uint8(1)
    
#             fx = np.asarray(
#                 self._update_network_synchronously_unchecked(x),
#                 dtype=np.uint8,
#             )
#             fy = np.asarray(
#                 self._update_network_synchronously_unchecked(y),
#                 dtype=np.uint8,
#             )
    
#             total_dist += float(np.sum(fx != fy))
    
#         return float(total_dist / float(n_simulations))

# # ===================== #
# #   Modular BoolForge   #
# # ===================== #

#     # def get_attractors_synchronous_exact_non_autonomous(self,
#     #     non_periodic_component : Sequence[Sequence[int]],
#     #     periodic_component : Sequence[Sequence[int]]) -> dict:
#     #     """
#     #     Compute all attractors and basin sizes under synchronous updating
#     #     for a Boolean network driven by a non-autonomous input sequence.
        
#     #     The input is split into a non-periodic component (applied once)
#     #     followed by a periodic component (repeated indefinitely). The
#     #     non-periodic component is first evaluated to determine a set of
#     #     initial states, which are then used to compute attractors under
#     #     the periodic component.
        
#     #     Parameters
#     #     ----------
#     #     non_periodic_component : sequence of sequence of int
#     #         External input values applied before the periodic regime.
#     #         Each inner sequence corresponds to one identity node and
#     #         contains binary values (0 or 1) over time.
        
#     #     periodic_component : sequence of sequence of int
#     #         External input values defining the periodic regime.
#     #         Each inner sequence corresponds to one identity node and
#     #         contains binary values (0 or 1) forming a repeating pattern.
        
#     #     Returns
#     #     -------
#     #     result : dict
#     #         Dictionary with the following keys:
        
#     #         - Attractors : list
#     #             List of attractors. Each attractor is a list of pairs
#     #             (external_input_decimal, state_decimal) forming a cycle.
        
#     #         - NumberOfAttractors : int
#     #             Total number of unique attractors.
        
#     #         - BasinSizes : list of int
#     #             Number of initial states converging to each attractor.
        
#     #         - AttractorDict : dict
#     #             Mapping from (external_input_decimal, state_decimal)
#     #             to attractor index.
        
#     #         - STG : dict
#     #             State transition graph mapping
#     #             (external_input_decimal, state_decimal) to the next pair.
        
#     #         - InitialStatesPeriodic : list of int
#     #             Initial state values (decimal) after applying the
#     #             non-periodic component.
        
#     #         - FormattedAttractors : list
#     #             Attractors represented as binary vectors, where the
#     #             external input bits and state bits are concatenated.
#     #     """
#     #     # Convert components into single argument? tuple|list|arr, str, etc.?
#     #     N = self.N - len(self.get_identity_nodes(False))
#     #     if len(non_periodic_component) > 0:
#     #         initial_states = set() # stores initial states for periodic computation
#     #         len_np_comp = len(non_periodic_component)
#     #         max_len_pattern = max(list(zip(map(len, non_periodic_component))))[0]
#     #         fixed_source_networks = {}
#     #         for i in range(2 ** N):
#     #             fxvec = utils.dec2bin(i, N) # initialize binary vector
#     #             for iii in range(max_len_pattern):
#     #                 values = [ non_periodic_component[j][iii] for j in range(len_np_comp) ]
#     #                 values_decimal = utils.bin2dec(values)
#     #                 if values_decimal in fixed_source_networks:
#     #                     fixed_source_network = fixed_source_networks[values_decimal]
#     #                 else:
#     #                     fixed_source_network = self.get_network_with_fixed_identity_nodes(values)
#     #                     fixed_source_networks[values_decimal] = fixed_source_network
#     #                 fxvec = fixed_source_network.update_network_synchronously(fxvec)
#     #             initial_states.add(utils.bin2dec(fxvec))
#     #         initial_states = list(initial_states)
#     #     else:
#     #         initial_states = list(range(2**N))
        
#     #     attr_computation = self.get_attractors_synchronous_exact_with_external_inputs(periodic_component, initial_states)
        
#     #     bvec_attractors = []
#     #     len_pattern = len(periodic_component)
#     #     for attr in attr_computation["Attractors"]:
#     #         bvec_attractors.append([])
#     #         for decimal_external, decimal_module in attr:
#     #             if len_pattern > 0:
#     #                 bvec = utils.dec2bin(decimal_external, len_pattern)
#     #             else:
#     #                 bvec = []
#     #             bvec.extend(utils.dec2bin(decimal_module, N))
#     #             bvec_attractors[-1].append(bvec)
        
#     #     attr_computation.update({"InitialStatesPeriodic":initial_states,"FormattedAttractors":bvec_attractors})
#     #     return attr_computation
    
#     # def get_attractors_synchronous_exact_with_external_inputs(self,
#     #     input_patterns : Sequence[Sequence[int]],
#     #     starting_states : [Sequence[int], None] = None) -> dict:
#     #     """
#     #     Compute all attractors and basin sizes under synchronous updating
#     #     for a Boolean network with periodic external inputs.
        
#     #     The external inputs are treated as a periodic sequence. The state
#     #     transition graph is constructed over the combined space of
#     #     (network state, input phase), and attractors are detected exactly.
        
#     #     Parameters
#     #     ----------
#     #     input_patterns : sequence of sequence of int
#     #         Periodic external input patterns. Each inner sequence
#     #         corresponds to one identity node and contains binary
#     #         values (0 or 1).
        
#     #     starting_states : sequence of int, optional
#     #         Optional list of initial network states in decimal form.
#     #         If None, all possible states are used.
        
#     #     Returns
#     #     -------
#     #     result : dict
#     #         Dictionary with the following keys:
        
#     #         - Attractors : list
#     #             List of attractors. Each attractor is a list of pairs
#     #             (external_input_decimal, state_decimal) forming a cycle.
        
#     #         - NumberOfAttractors : int
#     #             Total number of unique attractors.
        
#     #         - BasinSizes : list of int
#     #             Number of initial states converging to each attractor.
        
#     #         - AttractorDict : dict
#     #             Mapping from (external_input_decimal, state_decimal)
#     #             to attractor index.
        
#     #         - STG : dict
#     #             State transition graph mapping
#     #             (external_input_decimal, state_decimal) to the next pair.
#     #     """
#     #     N = self.N - len(self.get_identity_nodes(False))
        
#     #     if starting_states is None:
#     #         starting_states = list(range(2**N))
        
#     #     len_patterns = len(input_patterns)
#     #     lcm = math.lcm(*list(map(len, input_patterns)))
#     #     periodic_pattern_of_external_inputs = np.zeros((lcm, len_patterns), int)
#     #     for i, pattern in enumerate(input_patterns):
#     #         for j in range(int(lcm / len(pattern))):
#     #             periodic_pattern_of_external_inputs[len(pattern)*j:len(pattern)*(j+1),i] = pattern
#     #     n_initial_values = len(periodic_pattern_of_external_inputs)
        
#     #     fixed_source_networks = []
#     #     for input_values in periodic_pattern_of_external_inputs:
#     #         fixed_source_networks.append(self.get_network_with_fixed_identity_nodes(input_values))
        
#     #     lstt = utils.get_left_side_of_truth_table(N)
#     #     po2 = np.array([2**i for i in range(N)])[::-1]
        
#     #     dictF_fixed_source = []
        
#     #     for iii in range(n_initial_values):
#     #         state_space = np.zeros((2**N, N), dtype=int)
#     #         for i in range(N):
#     #             for j, x in enumerate(itertools.product([0, 1], repeat=fixed_source_networks[iii].indegrees[i])):
#     #                 if fixed_source_networks[iii].F[i][j]==1:
#     #                     # For rows in left_side_of_truth_table where the columns I[i] equal x, set state_space accordingly.
#     #                     state_space[np.all(lstt[:, fixed_source_networks[iii].I[i]] == np.array(x), axis=1), i] = 1
#     #         dictF_fixed_source.append(dict(zip(list(range(2**N)), np.dot(state_space, po2))))
        
#     #     attractors = []
#     #     basin_sizes = []
#     #     attractor_dict = dict()
#     #     stg = dict()
#     #     for iii_start in range(lcm):
#     #         for xdec in starting_states:
#     #             iii = iii_start
#     #             queue = [xdec]
#     #             while True:
#     #                 fxdec = dictF_fixed_source[iii % n_initial_values][xdec]
#     #                 stg.update({(int(utils.bin2dec(periodic_pattern_of_external_inputs[iii % n_initial_values])),int(xdec)):(int(utils.bin2dec(periodic_pattern_of_external_inputs[(iii + 1) % n_initial_values])),int(fxdec))})
#     #                 iii += 1
#     #                 try:
#     #                     index_attr = attractor_dict[(iii % n_initial_values,fxdec)]
#     #                     basin_sizes[index_attr] += 1
#     #                     attractor_dict.update(list(zip(zip(np.arange(iii_start,len(queue)+iii_start)%n_initial_values,queue), [index_attr] * len(queue))))
#     #                     break
#     #                 except KeyError:
#     #                     try: 
#     #                         index = queue[-n_initial_values::-n_initial_values].index(fxdec)
#     #                         dummy = np.arange(iii_start,len(queue)+iii_start)%n_initial_values
#     #                         #print(iii_start,j,list(zip(dummy[-n_initial_values*(index+1):],queue[-n_initial_values*(index+1):])))
#     #                         attractor_dict.update(list(zip(zip(dummy,queue), [len(attractors)] * len(queue))))
#     #                         attractors.append(list(zip(dummy[-n_initial_values*(index+1):],queue[-n_initial_values*(index+1):])))
#     #                         basin_sizes.append(1)
#     #                         break
#     #                     except ValueError:
#     #                         pass
#     #                 queue.append(fxdec)
#     #                 xdec = fxdec
        
#     #     attrs = []
#     #     attr_dict = {}
#     #     for key in attractor_dict.keys():
#     #         attr_dict[(int(utils.bin2dec(periodic_pattern_of_external_inputs[key[0]])), int(key[1]))] = int(attractor_dict[key])
#     #     for attr in attractors:
#     #         formatted_attr = []
#     #         for state in attr:
#     #             formatted_attr.append((int(utils.bin2dec(periodic_pattern_of_external_inputs[state[0]])), int(state[1])))
#     #         attrs.append(formatted_attr)
        
#     #     return { "Attractors":attrs, 
#     #             "NumberOfAttractors":len(attrs),
#     #             "BasinSizes":basin_sizes, 
#     #             "AttractorDict":attr_dict,
#     #             "STG":stg }


#     def _compute_post_transient_states(
#         self,
#         transient_input_sequence
#     ):
#         """
#         Apply the transient input sequence once to every state
#         and return the resulting unique states.
#         """
#         N_identity_nodes = len(self.get_identity_nodes(False))
#         N_regulated_nodes = self.N - N_identity_nodes
    
#         if not transient_input_sequence:
#             return list(range(2**N_regulated_nodes))
    
#         max_len = max(len(seq) for seq in transient_input_sequence)
#         num_inputs = len(transient_input_sequence)
    
#         fixed_network_cache = {}
#         resulting_states = set()
    
#         for state in range(2**N_regulated_nodes):
    
#             vec = utils.dec2bin(state, N_regulated_nodes)
    
#             for t in range(max_len):
#                 values = [
#                     transient_input_sequence[i][t]
#                     for i in range(num_inputs)
#                 ]
    
#                 values_dec = utils.bin2dec(values)
    
#                 if values_dec not in fixed_network_cache:
#                     fixed_network_cache[values_dec] = self.get_network_with_fixed_identity_nodes(values)
    
#                 vec = fixed_network_cache[values_dec].update_network_synchronously(vec)
    
#             resulting_states.add(utils.bin2dec(vec))
    
#         return list(resulting_states)


#     def _format_attractors_as_binary(
#         self,
#         attractors
#     ):
#         """
#         Convert attractors from (external_decimal, state_decimal)
#         to concatenated binary vectors.
#         """
    
#         N_identity_nodes = len(self.get_identity_nodes(False))
#         N_regulated_nodes = self.N - N_identity_nodes
    
#         formatted = []
    
#         for attr in attractors:
#             formatted_attr = []
    
#             for external_dec, state_dec in attr:
    
#                 if N_identity_nodes > 0:
#                     ext_bits = utils.dec2bin(external_dec, N_identity_nodes)
#                 else:
#                     ext_bits = []
    
#                 state_bits = utils.dec2bin(state_dec, N_regulated_nodes)
    
#                 formatted_attr.append(ext_bits + state_bits)
    
#             formatted.append(formatted_attr)
    
#         return formatted


#     def get_attractors_synchronous_exact_non_autonomous(self,
#         transient_input_sequence : Sequence[Sequence[int]],
#         periodic_input_sequence : Sequence[Sequence[int]]) -> dict:
#         """
#         Compute all attractors and basin sizes under synchronous updating
#         for a Boolean network driven by an external (non-autonomous) input
#         sequence.
        
#         The external input sequence is split into two parts:
        
#         1. A transient input sequence, applied once.
#         2. A periodic input sequence, repeated indefinitely.
        
#         First, the transient input sequence is applied to every network state
#         to determine the set of states that serve as initial conditions for
#         the periodic regime. Then, attractors and basin sizes are computed
#         exactly under the periodic input sequence.
        
#         Parameters
#         ----------
#         transient_input_sequence : sequence of sequence of int
#             External input values applied during the transient phase.
#             Each inner sequence corresponds to one external (identity)
#             node and contains binary values (0 or 1) indexed by time.
#             All sequences must have the same length.
        
#         periodic_input_sequence : sequence of sequence of int
#             External input values defining the periodic regime.
#             Each inner sequence corresponds to one external (identity)
#             node and contains binary values (0 or 1) forming a repeating
#             pattern. Each sequence may have different length; the
#             periodic regime is determined by their least common multiple.
        
#         Returns
#         -------
#         result : dict
#             Dictionary with the following keys:
        
#             - Attractors : list
#                 List of attractors. Each attractor is represented as a list
#                 of pairs (external_input_decimal, state_decimal) describing
#                 one full cycle in the augmented (input, state) space.
        
#             - NumberOfAttractors : int
#                 Total number of unique attractors.
        
#             - BasinSizes : list of int
#                 Number of initial states (after the transient phase) that
#                 converge to each attractor.
        
#             - AttractorDict : dict
#                 Mapping from (external_input_decimal, state_decimal) to
#                 attractor index.
        
#             - STG : dict
#                 State transition graph of the periodic regime, mapping
#                 (external_input_decimal, state_decimal) to its successor.
        
#             - InitialStatesPeriodic : list of int
#                 Set of network states (decimal representation) obtained
#                 after applying the transient input sequence.
        
#             - FormattedAttractors : list
#                 Attractors represented as binary vectors obtained by
#                 concatenating the external input bits and the network
#                 state bits at each point in the cycle.
#         """
    
#         #Apply transient block
#         initial_states = self._compute_post_transient_states(
#             transient_input_sequence
#         )
    
#         #Compute periodic attractors
#         result = self.get_attractors_synchronous_exact_with_external_inputs(
#             periodic_input_sequence,
#             initial_states
#         )
    
#         #Format attractors
#         result["InitialStatesPeriodic"] = initial_states
#         result["FormattedAttractors"] = self._format_attractors_as_binary(
#                                             result["Attractors"]
#                                         )
    
#         return result

#     def _build_periodic_input_matrix(self, periodic_input_sequence):
#         lengths = [len(p) for p in periodic_input_sequence]
#         lcm = math.lcm(*lengths)
    
#         periodic_inputs = np.zeros((lcm, len(periodic_input_sequence)), dtype=int)
    
#         for i, pattern in enumerate(periodic_input_sequence):
#             reps = lcm // len(pattern)
#             periodic_inputs[:, i] = pattern * reps
    
#         return periodic_inputs, lcm
    
#     def _build_phase_transition_maps(self, periodic_inputs):
#         N_regulated_nodes = self.N - len(self.get_identity_nodes(False))
    
#         transition_maps = []
    
#         for phase_values in periodic_inputs:
#             fixed_net = self.get_network_with_fixed_identity_nodes(phase_values)
    
#             phase_map = {}
#             for state in range(2**N_regulated_nodes):
#                 next_state = utils.bin2dec(
#                     fixed_net.update_network_synchronously(
#                         utils.dec2bin(state, N_regulated_nodes)
#                     )
#                 )
#                 phase_map[state] = next_state
    
#             transition_maps.append(phase_map)
    
#         return transition_maps
    
#     def _explore_phase_state_space(self, transition_maps, periodic_inputs, starting_states_dec):
#         lcm = len(periodic_inputs)
    
#         attractors = []
#         basin_sizes = []
#         attractor_dict = {}
#         stg = {}
        
#         for phase in range(lcm):
#             for state in starting_states_dec:
#                 path = []
#                 visited_local = {}
    
#                 while True:
    
#                     key = (phase, state)
    
#                     if key in attractor_dict:
#                         idx = attractor_dict[key]
#                         basin_sizes[idx] += 1
#                         break
    
#                     if key in visited_local:
#                         cycle_start = visited_local[key]
#                         cycle = path[cycle_start:]
#                         idx = len(attractors)
    
#                         attractors.append(cycle)
#                         basin_sizes.append(1)
    
#                         for k in cycle:
#                             attractor_dict[k] = idx
#                         break
    
#                     visited_local[key] = len(path)
#                     path.append(key)
    
#                     next_state = transition_maps[phase][state]
#                     next_phase = (phase + 1) % lcm
    
#                     stg[(utils.bin2dec(periodic_inputs[phase]), state)] = (
#                         utils.bin2dec(periodic_inputs[next_phase]),
#                         next_state
#                     )
    
#                     phase = next_phase
#                     state = next_state
    
#         return attractors, basin_sizes, attractor_dict, stg    
    
    
#     def get_attractors_synchronous_exact_with_external_inputs(
#         self,
#         periodic_input_sequence : Sequence[Sequence[int]],
#         starting_states_dec : [Sequence[int], None] = None) -> dict:
#         """
#         Compute all attractors and basin sizes under synchronous updating
#         for a Boolean network with periodic external inputs.
        
#         The external inputs are treated as a periodic sequence. The state
#         transition graph is constructed over the combined space of
#         (network state, input phase), and attractors are detected exactly.
        
#         Parameters
#         ----------
#         periodic_input_sequence : sequence of sequence of int
#             External input values defining the periodic regime.
#             Each inner sequence corresponds to one external (identity)
#             node and contains binary values (0 or 1) forming a repeating
#             pattern. Each sequence may have different length; the
#             periodic regime is determined by their least common multiple.
        
#         starting_states_dec : sequence of int, optional
#             Optional list of initial network states in decimal representation.
#             If None, all possible states are used.
        
#         Returns
#         -------
#         result : dict
#             Dictionary with the following keys:
        
#             - Attractors : list
#                 List of attractors. Each attractor is a list of pairs
#                 (external_input_decimal, state_decimal) forming a cycle.
        
#             - NumberOfAttractors : int
#                 Total number of unique attractors.
        
#             - BasinSizes : list of int
#                 Number of initial states converging to each attractor.
        
#             - AttractorDict : dict
#                 Mapping from (external_input_decimal, state_decimal)
#                 to attractor index.
        
#             - STG : dict
#                 State transition graph mapping
#                 (external_input_decimal, state_decimal) to the next pair.
#         """    
#         N_regulated_nodes = self.N - len(self.get_identity_nodes(False))
    
#         if starting_states_dec is None:
#             starting_states_dec = list(range(2**N_regulated_nodes))
    
#         periodic_inputs, lcm = self._build_periodic_input_matrix(periodic_input_sequence)
    
#         transition_maps = self._build_phase_transition_maps(periodic_inputs)
    
#         attractors, basin_sizes, attractor_dict, stg = \
#             self._explore_phase_state_space(
#                 transition_maps,
#                 periodic_inputs,
#                 starting_states_dec
#             )
    
#         formatted_attractors = []
#         for attr in attractors:
#             formatted_attractors.append([
#                 (utils.bin2dec(periodic_inputs[phase]), state)
#                 for phase, state in attr
#             ])
    
#         return {
#             "Attractors": formatted_attractors,
#             "NumberOfAttractors": len(formatted_attractors),
#             "BasinSizes": basin_sizes,
#             "AttractorDict": attractor_dict,
#             "STG": stg
#         }

#     def _get_fnet_(self,values,fixed_network_cache):
#         values_dec = utils.bin2dec(values)
#         if values_dec in fixed_network_cache:
#             fixed_network = fixed_network_cache[values_dec]
#         else:
#             fixed_network = self.get_network_with_fixed_identity_nodes(values)
#             fixed_network_cache[values_dec] = fixed_network
#         return fixed_network

#     def _calculate_trajectory(
#             self,
#             starting_state_dec,
#             transient_input_sequence,
#             periodic_input_sequence,
#             N_regulated_nodes,
#             fixed_network_cache
#             ):
#         trajectory = [starting_state_dec]
#         latest_state = starting_state_dec
        
#         # Compute the non-periodic component of the trajectory.
#         len_np = len(transient_input_sequence)
#         max_len_pattern = max(list(zip(map(len, transient_input_sequence))))[0]
#         for idx in range(max_len_pattern):
#             vals = [ transient_input_sequence[node][idx] for node in range(len_np) ]
#             fixed_network = self._get_fnet_(vals,fixed_network_cache)
#             latest_state = utils.bin2dec(
#                 fixed_network.update_network_synchronously(
#                     utils.dec2bin(latest_state, N_regulated_nodes)
#                 )
#             )
#             trajectory.append(latest_state)
            
#         # Compute the periodic component of the trajectory.
#         len_p = len(periodic_input_sequence)
#         lcm = math.lcm(*list(map(len, periodic_input_sequence)))
#         idx_p = 0
#         cycle_len = -1
        
#         seen = {}  # (state, phase) -> index in traj_cyclic
#         traj_cyclic = []
#         idx_p = 0
        
#         while True:
#             phase = idx_p % lcm
#             key = (latest_state, phase)
        
#             if key in seen:
#                 # We found the cycle start
#                 cycle_start = seen[key]
#                 cycle_len = len(traj_cyclic) - cycle_start
#                 break
        
#             seen[key] = len(traj_cyclic)
        
#             vals = [
#                 periodic_input_sequence[node][phase]
#                 for node in range(len_p)
#             ]
#             fixed_network = self._get_fnet_(vals,fixed_network_cache)
        
#             latest_state = utils.bin2dec(
#                 fixed_network.update_network_synchronously(
#                     utils.dec2bin(latest_state, N_regulated_nodes)
#                 )
#             )
        
#             traj_cyclic.append(latest_state)
#             idx_p += 1
#         trajectory.extend(traj_cyclic[:cycle_start + cycle_len])
#         #print(trajectory, traj_cyclic, traj_cyclic[:cycle_start + cycle_len])
        
#         # Compress the trajectory's representation to be minimal.
#         # That is, only the transient input sequence and a single
#         # cycle of the periodic input sequence.
#         len_traj = len(trajectory)
#         best_trajectory = []
#         best_cycle_len = -1
#         best_length = math.inf
#         for s in range(len_traj):
#             for p in range(1, min(cycle_len, len_traj - s) + 1):
#                 proposed_period = trajectory[s : s + p]
#                 good_proposal = True
#                 for i in range(s, len_traj):
#                     if trajectory[i] != proposed_period[(i - s) % p]:
#                         good_proposal = False
#                         break
#                 if not good_proposal:
#                     continue
                
#                 len_proposal = s + p
#                 if len_proposal < best_length:
#                     best_length = len_proposal
#                     best_trajectory = trajectory[:s] + proposed_period
#                     best_cycle_len = p
#         #print(best_trajectory, best_cycle_len, "\n")
        
#         # Return the compressed trajectory array and the length of the
#         # periodic component.
#         # Note that the periodic_input_sequence will ALWAYS be the last
#         # cycle_len values in the array. The periodic_input_sequence
#         # also correspond with the attractors of the network.
#         return best_trajectory, best_cycle_len
    

    
    
#     def get_trajectories(
#         self,
#         transient_input_sequence,
#         periodic_input_sequence,
#         merge_trajectories=True,
#         starting_states_dec=None
#     ):
#         """
#         Compute synchronous trajectories of the Boolean network given a 
#         non-autonomous external input sequence.
        
#         The external input is split into two phases:
        
#         1. A transient input sequence, applied once.
#         2. A periodic input sequence, repeated indefinitely.
        
#         For each specified initial state, the transient input sequence is
#         applied first. The resulting state then evolves under the periodic
#         input sequence. Periodicity is detected in the augmented space
#         (state, input phase), ensuring that cycles are identified correctly
#         even when the same network state appears at different input phases.
        
#         Each trajectory is returned in minimal form, consisting of:
        
#         - A non-periodic prefix (possibly empty),
#         - Followed by a single instance of the periodic cycle.
        
#         Parameters
#         ----------
#         transient_input_sequence : sequence of sequence of int
#             External input values applied during the transient phase.
#             Each inner sequence corresponds to one external (identity)
#             node and contains binary values (0 or 1) indexed by time.
#             All sequences must have the same length.
        
#         periodic_input_sequence : sequence of sequence of int
#             External input values defining the periodic regime.
#             Each inner sequence corresponds to one external (identity)
#             node and contains binary values (0 or 1) forming a repeating
#             pattern. The effective period is the least common multiple
#             of the individual sequence lengths.
        
#         merge_trajectories : bool, optional (default=True)
#             If True, trajectories are merged into a directed graph
#             representing the non-autonomous state space (with consistent
#             merging across trajectories). If False, individual trajectories
#             are returned.
        
#         starting_states_dec : sequence of int or None, optional
#             Decimal representations of initial network states.
#             If None, all states in the full state space (size 2^N) are used.
        
#         Returns
#         -------
#         result : networkx.DiGraph or list
#             If merge_trajectories is True:
#                 A directed graph representing the merged trajectories,
#                 i.e., the state space of the non-autonomous system.
        
#             If merge_trajectories is False:
#                 A list of tuples (trajectory, cycle_length), where
        
#                 - trajectory : list of int
#                     Decimal representations of states, containing the
#                     non-periodic prefix followed by exactly one full
#                     cycle of the periodic component.
        
#                 - cycle_length : int
#                     Length of the periodic component (the last
#                     cycle_length entries of trajectory).
#         """        
#         N_identity_nodes = len(self.get_identity_nodes(False))
#         N_regulated_nodes = self.N - N_identity_nodes
    
#         # Validation
#         assert len(transient_input_sequence) == len(periodic_input_sequence)
#         assert len(transient_input_sequence) == N_identity_nodes
#         assert all(len(seq) > 0 for seq in periodic_input_sequence)
    
#         if starting_states_dec is None:
#             starting_states_dec = list(range(2 ** N_regulated_nodes))
#         else:
#             starting_states_dec = list(set(starting_states_dec))
    
#         fixed_network_cache = {}
#         trajectories = []
    
#         for state in starting_states_dec:
#             trajectory, cycle_len = self._calculate_trajectory(
#                 state,
#                 transient_input_sequence,
#                 periodic_input_sequence,
#                 N_regulated_nodes,
#                 fixed_network_cache
#             )        
#             trajectories.append((trajectory, cycle_len))
    
#         if merge_trajectories:
#             return utils.compress_trajectories(trajectories, N_regulated_nodes)
    
#         return trajectories
    