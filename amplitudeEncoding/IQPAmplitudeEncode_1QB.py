'''
We alter the IQP ansatz to include Amplitude Encoding onto one qubit. 
'''

from __future__ import annotations

# __all__ = ['CircuitAnsatz',
#            'IQPAnsatz',
#            'Sim14Ansatz',
#            'Sim15Ansatz',
#            'StronglyEntanglingAnsatz']

from abc import abstractmethod
from collections.abc import Callable, Mapping
from itertools import cycle
from typing import Type

import numpy as np
from sympy import Symbol, symbols

from lambeq.ansatz import BaseAnsatz
from lambeq.backend.grammar import Box, Diagram, Functor, Ty
from lambeq.backend.quantum import (
    Bra,
    CRz,
    Diagram as Circuit,
    Discard,
    H,
    Id,
    Ket,
    quantum,
    qubit,
    Rotation,
    Rx, Ry, Rz
)

from lambeq import CircuitAnsatz

computational_basis = Id(qubit)

from lambeq import AtomicType

from lambeq import BobcatParser

N = AtomicType.NOUN
S = AtomicType.SENTENCE
adj = N @ N.l

parser = BobcatParser(verbose='text')


'''
General Ansatz class that we changed for our purpose
'''

class CircuitAnsatz(BaseAnsatz):
    """Base class for circuit ansatz."""

    def __init__(self,
                 ob_map: Mapping[Ty, int],
                 n_layers: int,
                 n_single_qubit_params: int,
                 circuit: Callable[[int, np.ndarray], Circuit],
                 nounParams: np.ndarray, 
                 discard: bool = False,
                 single_qubit_rotations: list[Type[Rotation]] | None = None,
                 postselection_basis: Circuit = computational_basis) -> None:
        """Instantiate a circuit ansatz.

        Parameters
        ----------
        ob_map : dict
            A mapping from :py:class:`lambeq.backend.grammar.Ty` to
            the number of qubits it uses in a circuit.
        n_layers : int
            The number of layers used by the ansatz.
        n_single_qubit_params : int
            The number of single qubit rotations used by the ansatz.
        circuit : callable
            Circuit generator used by the ansatz. This is a function
            (or a class constructor) that takes a number of qubits and
            a numpy array of parameters, and returns the ansatz of that
            size, with parameterised boxes.
        discard : bool, default: False
            Discard open wires instead of post-selecting.
        postselection_basis: Circuit, default: Id(qubit)
            Basis to post-select in, by default the computational basis.
        single_qubit_rotations: list of Circuit, optional
            The rotations to be used for a single qubit. When only a
            single qubit is present, the ansatz defaults to applying a
            series of rotations in a cycle, determined by this parameter
            and `n_single_qubit_params`.

        """
        self.ob_map = {src: qubit ** ty if isinstance(ty, int) else ty
                       for src, ty in ob_map.items()}
        self.n_layers = n_layers
        self.n_single_qubit_params = n_single_qubit_params
        self.circuit = circuit
        self.nounParams = nounParams
        if len(nounParams) != n_single_qubit_params: 
            print('The number of noun parameters does not fit the number of single qubit operators. ')
        self.discard = discard
        self.postselection_basis = postselection_basis
        self.single_qubit_rotations = single_qubit_rotations or []

        self.functor = Functor(target_category=quantum,
                               ob=self._ob,
                               ar=self._ar)


    def __call__(self, diagram: Diagram) -> Circuit:
        """Convert a lambeq diagram into a lambeq circuit."""
        return self.functor(diagram)  # type: ignore[return-value]


    def ob_size(self, pg_type: Ty) -> int:
        """Calculate the number of qubits used for a given type."""
        return sum(map(len, map(self.functor, pg_type)))


    @abstractmethod
    def params_shape(self, n_qubits: int) -> tuple[int, ...]:
        """Calculate the shape of the parameters required."""


    def _ob(self, _: Functor, ty: Ty) -> Ty:
        return self.ob_map[ty]

    def _ar(self, _: Functor, box: Box) -> Circuit:
        label = self._summarise_box(box)
        dom, cod = self.ob_size(box.dom), self.ob_size(box.cod)

        n_qubits = max(dom, cod)
        if n_qubits == 0:
            circuit = Id()

        # NOUNS AND OTHER OPERATIONS THAT ACT ON ONE QUBIT ONLY
        elif n_qubits == 1:

            syms = symbols(f'{label}_0:{self.n_single_qubit_params}',
                           cls=Symbol)

            circuit = Id(qubit)
            # for rot, sym in zip(cycle(self.single_qubit_rotations), syms):
            #     print('rot', rot)
            #     circuit >>= rot(self.nounParams[0])

            # cycle ensures that we use all three parameters that we are given, because there is only two parameters specified! Rx and Rz, but we want Rx Rz Rx. 
            # for rot, nounParam in zip(cycle(self.single_qubit_rotations), self.nounParams): 
            #     circuit >>= rot(nounParam)

            if 'sauce' in str(syms[0]) or 'meal' in str(syms[0]) or 'dinner' in str(syms[0]) or 'chef' in str(syms[0]):
                circuit >>= Ry(0)
                print('we choose the paramter 0 for the symbol: ', str(syms[0]))
            elif 'program' in str(syms[0]) or 'application' in str(syms[0]) or 'software' in str(syms[0]) or 'programmer' in str(syms[0]):
                if '†' in str(syms[0]):
                    circuit >>= Ry(-np.pi)
                else: 
                    circuit >>= Ry(np.pi)
                print('we choose the paramter pi for the symbol: ', str(syms[0]))
            elif 'woman' in str(syms[0]) or 'man' in str(syms[0]) or 'person' in str(syms[0]):
                print('we choose the paramter pi/2 for the symbol: ', str(syms[0]))
                if '†' in str(syms[0]):
                    circuit >>= Ry(-np.pi/2)
                else: 
                    circuit >>= Ry(np.pi/2)
        else:
            params_shape = self.params_shape(n_qubits)

            # Verbs and other non-nouns

            syms = symbols(f'{label}_0:{np.prod(params_shape)}', cls=Symbol)

            # indicesToDelete = []
            # symsToChange = []
            # 
            # for sym_no in range(len(syms)):
            #     sym = syms[sym_no]
            #     print(sym)
            #     if '_n_' in str(sym): 
            #         # indicesToDelete.append(sym_no)
            #         symsToChange.append(sym) 
            # 
            # syms = [syms[i] for i in range(len(syms)) if i not in indicesToDelete]

            
            params: np.ndarray = np.array(syms).reshape(params_shape)
            circuit = self.circuit(n_qubits, params)

        if cod > dom:
            circuit = Id(dom) @ Ket(*[0]*(cod - dom)) >> circuit
        elif cod < dom:
            if self.discard:
                circuit >>= Id(cod) @ Id().tensor(
                    *[Discard() for _ in range(dom - cod)]
                )
            else:
                circuit >>= Id(cod).tensor(
                    *[self.postselection_basis] * (dom-cod))
                circuit >>= Id(cod) @ Bra(*[0]*(dom - cod))
        return circuit
    








class IQPAmplitudeEncode(CircuitAnsatz):
    """Instantaneous Quantum Polynomial ansatz.

    An IQP ansatz interleaves layers of Hadamard gates with diagonal
    unitaries. This class uses :py:obj:`n_layers-1` adjacent CRz gates
    to implement each diagonal unitary.

    Code adapted from DisCoPy.

    """

    def __init__(self,
                 ob_map: Mapping[Ty, int],
                 n_layers: int,
                 nounParams: np.ndarray,
                 n_single_qubit_params: int = 3,
                 discard: bool = False) -> None:
        """Instantiate an IQP ansatz.

        Parameters
        ----------
        ob_map : dict
            A mapping from :py:class:`lambeq.backend.grammar.Ty` to
            the number of qubits it uses in a circuit.
        n_layers : int
            The number of layers used by the ansatz.
        n_single_qubit_params : int, default: 3
            The number of single qubit rotations used by the ansatz.
        discard : bool, default: False
            Discard open wires instead of post-selecting.

        """
        super().__init__(ob_map,
                         n_layers,
                         n_single_qubit_params,
                         self.circuit,
                         nounParams, 
                         discard,
                         [Rx, Rz])


    def params_shape(self, n_qubits: int) -> tuple[int, ...]:
        return (self.n_layers, n_qubits - 1)


    def circuit(self, n_qubits: int, params: np.ndarray) -> Circuit:
        
        if n_qubits == 1:
            # circuit = Rx(params[0]) >> Rz(params[1]) >> Rx(params[2])
            print('test')
        else:
            circuit = Id(n_qubits)
            hadamards = Id().tensor(*(n_qubits * [H]))
            for thetas in params:
                rotations = Id(n_qubits).then(*(
                    Id(i) @ CRz(thetas[i]) @ Id(n_qubits - 2 - i)
                    for i in range(n_qubits - 1)))
                circuit >>= hadamards >> rotations
            if self.n_layers > 0:  # Final layer of Hadamards
                circuit >>= hadamards

        return circuit  # type: ignore[return-value]


