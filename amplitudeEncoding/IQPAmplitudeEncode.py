'''
In this file, we capture the rewritten ansatz to amplitude encode nouns in a circuit. 
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
    Rx, Ry, Rz, 
    CRy, 
    X,
    SWAP
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
Functions that we need for the Ansatz 
'''

def noun_ansatz(paramString: str, parameterDict, circuitt):
    '''
    this function returns a 2-qubit circuit encoding the four dimensional vector arr 
    '''
    word, dagger = GetWordFromParam(paramString)
    print('the word is:  ', word)
    print('is it daggered? ', dagger)

    arr = parameterDict[word]

    print('and the parameters are: ', arr)

    if len(arr) != 4: 
        print('Error, we need a four dimensional vector')
    
    if arr == [0,0,0,0]:
        return circuitt

    for i in range(len(arr)): 
        if arr[i] == 0: 
            # avoid nan, I suppose? 
            arr[i] = 1e-10

    # normalize vector 
    arr = arr / np.linalg.norm(arr)

    print('vector to amplitude encode: ', arr)


    a1 = np.linalg.norm(arr[0:2])
    a2 = np.linalg.norm(arr[2:])
    phi1 = np.arccos(a1)/np.pi

    # fix issues with rotations
    rot1 = arr[0:2]/a1
    phi2_cos = np.arccos(rot1[0])/np.pi
    phi2_sin = np.arcsin(rot1[1])/np.pi
    if not np.sign(phi2_cos) == np.sign(phi2_sin):
        phi2_cos *= -1
    rot2 = arr[2: ]/a2
    phi3_cos = np.arccos(rot2[0])/np.pi
    phi3_sin = np.arcsin(rot2[1])/np.pi
    if not np.sign(phi3_cos) == np.sign(phi3_sin):
        phi3_cos *= -1

    # print('is this executed succesfully?? The parameters are: ', phi1)

    partToAdd = Ry(phi1) @ Id(1) >> CRy(phi3_cos) >> X @ Id(1) >> CRy(phi2_cos) >> X @ Id(1)

    # if the word is upside down, we flip the whole circuit 
    # print('debugdebugdebug, the thing we are going to add: ')
    # partToAdd.draw()
    
    if dagger:
        partToAddNew = SWAP >> partToAdd
        # print('and the corresponding dagger circuit: ')
        # partToAddNew.dagger().draw()
        circuitt >>= partToAddNew.dagger()

    # if the word is not upside down, we just append the circuit as is 
    else: 
        circuitt >>= partToAdd

    # making sure the amplitude encoding works as intended 

    # circuit.draw()
    # amplitude = circuit.eval()
    # print('the amplitude is: ', amplitude)
    # probability = abs(amplitude) ** 2   

    # print('the probability is: ', probability)

    # return circuit >> Ry(phi1) @ circuit.id(1) >> CRy(phi3_cos) >> X @ circuit.id(1) >> CRy(phi2_cos) >> X @ circuit.id(1)
    return circuitt

def GetWordFromParam(param: str):
    word = param.split('_')[0] 

    if '†' in word: 
        dagger = True 
        word = word.replace('†', '')
    else: 
        dagger = False 
    
    return word, dagger

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
                 parameterDict: dict, 
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
        self.parameterDict = parameterDict
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

            circuit >>= Ry(self.nounParams[0])

            # remove first element from list! 
            print('REMOVEREMOVE JKJAJDFJ')
            self.nounParams = self.nounParams[1:]

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
    


'''
This is a child of CircuitAnsatz, executing the amplitude encoding
'''

class IQPAmplitudeEncode2QB(CircuitAnsatz):
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
                 parameterDict: dict, 
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
                         parameterDict, 
                         discard,
                         [Rx, Rz])


    def params_shape(self, n_qubits: int) -> tuple[int, ...]:
        return (self.n_layers, n_qubits - 1)


    def circuit(self, n_qubits: int, params: np.ndarray) -> Circuit:

        # when we amplitude encode the nouns, we don't want to draw the final Hadamard layer 
        addH = True
        
        if n_qubits == 1:
            # circuit = Rx(params[0]) >> Rz(params[1]) >> Rx(params[2])
            print('test')
        else:
            circuit = Id(n_qubits)
            hadamards = Id().tensor(*(n_qubits * [H]))
            print('params in circuitfunction:', params)
            for thetas in params:

                # print('first draw of circuit: ')
                # circuit.draw()
                print('thetas: ', thetas)
                if '_n_' in str(thetas[0]): 
                    print('we have to do somehting')
                    # continue 
                    circuit = noun_ansatz(str(thetas[0]), self.parameterDict, circuitt=circuit)

                    # remove the used noun parameters 
                    # self.nounParams = self.nounParams[4:]

                    addH = False
                    continue
                rotations = Id(n_qubits).then(*(
                    Id(i) @ CRz(thetas[i]) @ Id(n_qubits - 2 - i)
                    for i in range(n_qubits - 1)))
                circuit >>= hadamards >> rotations
                # print('after hadamard and rotations: ')
                # circuit.draw()

            if self.n_layers > 0 and addH == True:  # Final layer of Hadamards
                circuit >>= hadamards

        return circuit  # type: ignore[return-value]










