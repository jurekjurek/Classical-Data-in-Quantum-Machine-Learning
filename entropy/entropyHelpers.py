'''
Function to go through the test diagrams one by one 
and modify them accordingly. 
'''
import numpy as np
from qiskit.circuit import Parameter
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.converters import dag_to_circuit, circuit_to_dag

def SetUpDiags(circuit: QuantumCircuit, parameterDict: dict, amplitudeEncoded: bool, WordsToForget: list = ['woman', 'person', 'meal'], oneQubit = True, traceOutInsteadOfMeasure:bool = True, show: bool = True):
    '''
    Given a circuit as an Input, we remove the noun parameters belonging to the first noun in the sentence. 

    We then, depending on the shape of the circuit, replace the noun parameters by an entangled Bell state or we forget about it. 

    In the end, we replace the verb parameters by the parameters learned by the model. 

    The reason we are doing this is to introduce a model with real uncertainty. 

    wordsToForget: list of strings of words that we want to forget ('man', 'woman', 'person', ...)

    '''
    # remove noun parameters and set the rest of the parameters 

    # start by extracting circuit data 
    circuitData = circuit.data

    if show:     
        print('circuit that we start with: ')
        print(circuit)


    # first, figure out the qubit the sentence is encoded on 
    measuredQubits = []
    for i in range(len(circuitData)):
            if str(circuitData[i].operation.name) == 'measure':
                measuredQubits.append(circuitData[i].qubits[0].index)

    qubitList = list(range(circuit.num_qubits))
    for qubit in qubitList:
        if qubit not in measuredQubits: 
            sentenceQubit = qubit 
            break 


    print('we determined the sentencqubit to be: ', sentenceQubit)


    # list to store the nounParameters to be removed 
    indexList = []

    # boolean that tells us if the noun of interest is daggered or not 
    dagger = False 

    for i in range(len(circuitData)): 

        if 'programmer' in str(circuitData[i].operation.params) and 'program' in WordsToForget: 
            continue

        # check if this gate has a parameter of a word we want to forget about
        for word in WordsToForget:
            
            if word in str(circuitData[i].operation.params):
                # store the parameter in a list
                indexList.append(i)

                print('we forget the following word: ', word)

                nounQubit = circuitData[i].qubits[0].index

                if '†' in str(circuitData[i].operation.params):
                    dagger = True 
                    
                # break


    # FIRST: UPDATE THE PARAMETERS, THEN REMOVE THE GATES FROM THE CIRCUIT 

    # get parameters of the remaining circuit 
    params = circuit.parameters

    print('AMPLITUDE ENCODING IS: ', amplitudeEncoded)

    # set the parameters of the remaining words based on the models weights 
    for i in range(len(params)):
        # print(str(params[i]))
        # print(parameterDict[str(params[i])])
        print(params[i])
        if str(params[i]) in parameterDict:

            print('DEBUGDEBUGDEBUGDEBUG,  if we put in parameter ', parameterDict[str(params[i])], 'for the label: ', str(params[i]))
            circuit = circuit.assign_parameters({params[i]: (parameterDict[str(params[i])])})#/(2*np.pi))})

        # else: 
        #     '''
        #     assign the parameter ... 
        #     '''


        elif amplitudeEncoded: 

            print('AEAEAEAEAEAEAEA!!!')
            if 'man' in str(params[i]) or 'woman' in str(params[i]) or 'person' in str(params[i]):

                # if 1 in params, this means that it is the x gate in between the z gates 

                if '†' in str(params[i]) and '0' in str(params[i]): 
                    circuit = circuit.assign_parameters({params[i]: -np.pi/2})
                elif '1' in str(params[i]):
                    circuit = circuit.assign_parameters({params[i]: np.pi/2})
                elif '0' in str(params[i]): 
                    circuit = circuit.assign_parameters({params[i]: np.pi/2})
                else: 
                    circuit = circuit.assign_parameters({params[i]: 0})
            
            elif 'software' in str(params[i]) or 'application' in str(params[i]) or 'program' in str(params[i]) or 'programmer' in str(params[i]): 
                if '0' in str(params[i]):
                    circuit = circuit.assign_parameters({params[i]: np.pi})
                else: 
                    circuit = circuit.assign_parameters({params[i]: 0})

            elif 'meal' in str(params[i]) or 'sauce' in str(params[i]) or 'dinner' in str(params[i]) or 'chef' in str(params[i]): 
                circuit = circuit.assign_parameters({params[i]: 0})


        # probably unnecesary
        else: 
            print('TESTSETESTSETS')
            # if the model only learned parameters for the word that are not daggered
            if '†' in str(params[i]):
                if oneQubit == True: 
                    tempParam = str(params[i])
                    # reverse string 
                    tempParam = tempParam[::-1]

                    # remove elements 2,3,4,5
                    tempParam = tempParam[:2] + 'n__' + tempParam[5:]

                    tempParam = tempParam[::-1]
                    print(tempParam)             

                    circuit = circuit.assign_parameters({params[i]: (parameterDict[tempParam])})#/(2*np.pi))}) 
                if oneQubit == False:
                    tempParam = str(params[i])
                    # reverse string 
                    tempParam = tempParam[::-1]

                    # remove elements 2,3,4,5
                    tempParam = tempParam[:2] + 'n__' + tempParam[5:]

                    tempParam = tempParam[::-1]
                    print(tempParam)             

                    circuit = circuit.assign_parameters({params[i]: (parameterDict[tempParam])})#/(2*np.pi))}) 

    if show: 
        print('Circuit with new parameters: ')
        print(circuit)
    

    if len(indexList) != 0: 
        # remove the parameters of the words we want to forget
        circuit.data.pop(indexList[0])
        circuit.data.pop(indexList[1]-1)
        circuit.data.pop(indexList[2]-2)

    if show: 
        print('circuit with removed words to forget: ')
        print(circuit)




    # if the word we want to forget about is an effect, we simply forget about the qubit 
    if dagger: 
        # remove measurement 
        circuitData = circuit.data

        for i in range(len(circuitData)):

            if str(circuitData[i].operation.name) == 'measure': 
                if circuitData[i].qubits[0].index == nounQubit: 

                    # remove measurement
                    circuit.data.pop(i)
                    break
        
        # that's it 
        
    

    # if the word we want to forget about is a state, we introduce a bell state and clip it to the beginning of the circuit 
    if not dagger: 
        # introduce circuit with 5 qubits to be composed with the initial circuit 
        bellCircuit = QuantumCircuit(circuit.num_qubits + 1)
        # bell state on first two qubits 
        bellCircuit.h(0)

        # both noun and sentencequbit move up one position, since we add a qubit in position 0
        nounQubit = nounQubit + 1
        sentenceQubit = sentenceQubit + 1

        # entangle the artificial qubit with the qubit of interest
        bellCircuit.cx(0,nounQubit)

        # now, add one 'artificial' qubit to the initial circuit 
        q = QuantumRegister(1, 'q')
        circuit = circuit.reverse_bits()
        circuit.add_bits(q)
        circuit = circuit.reverse_bits()

        # now that the two circuits are equal in length, combine them 
        circuit = bellCircuit.compose(circuit) #, qubits=[0,1,2,3,4])

        # now, we have our circuit 






    if show: 
        print('circuit with discarded qubit: ')
        print(circuit)

    circuit.draw(output= 'mpl', filename='testtest.png')



    # remove the measurements and just trace over the corresponding later 
    if traceOutInsteadOfMeasure: 
        circuitData = circuit.data

        i = len(circuitData) - 1

        traceList = []
        # for i in range(len(circuitData)):
        #     if str(circuitData[i].operation.name) == 'measure':
        #         traceList.append(circuitData[i].qubits[0].index)

        while i >= 0: 
            if str(circuitData[i].operation.name) == 'measure': 
                print('removed measurement. ')

                # if traceList over all -2 qubits 
                # traceList.append(circuitData[i].qubits[0].index)
                del circuit.data[i]
            i -= 1

        # if tracelist over all -1 qubits 
        nQ = circuit.num_qubits
        # traceList = [x for x in range(nQ) if x != sentenceQubit]

        # if dagger, we did not add anything, so we encode the meaning of the sentence on the 4x4 density matrix on the noun and sentence wire 
        if dagger: 
            traceList = [x for x in range(nQ) if (x != sentenceQubit and x != nounQubit)]

        # if not dagger, we have a bell circuit and the open qubit wire is in qubit 0s position 
        elif not dagger: 
            # bell added 
            traceList = [x for x in range(1, nQ) if x != sentenceQubit]

        print('trace out instead of measure, circuit is: ')
        print(circuit)
        print('And traceList is: ', traceList)


    


    if len(indexList) == 0: 
        print('nothing to remove.')

        # if nothing was to remove, nothing is going to be added, so addedQubit: bool = False
        return circuit, False


    if dagger:
        qubitToTraceOut = [nounQubit]
    if not dagger: 
        qubitToTraceOut = [0]

    return circuit, not dagger, traceList, qubitToTraceOut


from lambeq import BobcatParser
from lambeq import RemoveCupsRewriter
from lambeq import IQPAnsatz

from pytket.extensions.qiskit import tk_to_qiskit





'''
Now, for the density matrices 
'''
from qiskit import QuantumCircuit
from qiskit.quantum_info import DensityMatrix, partial_trace
from qiskit.visualization import plot_state_city



def GetDensityMatrix(circuit: QuantumCircuit, traceList: list, addedQubit: bool, traceOutInsteadOfMeasure: bool, qubitToTraceOut):
    '''
    Given a QuantumCircuit representing a sentence, this function returns the density matrix of the *sentence* qubit 
    ''' 
    import qiskit.quantum_info as qi 
    # from qiskit.providers.aer import AerSimulator
    
    # dmSimulation = AerSimulator()

    # job = dmSimulation.run(circuit.decompose(reps = 1))



    # result = job.result()

    # # print(result)

    # data = result.data()

    # rho = data.get('density_matrix')

    rho = DensityMatrix(circuit)


    # print('CIRCUIT IN DENSITY MATRIX')
    # print(circuit)

    # differentiate between dagger and not dagger (one has one qubit less)
    # dagger
    # if circuit.num_qubits == 4: 
    #     sentenceDM = partial_trace(state=rho, qargs=[0, 2, 3])

    # # not dagger 
    # elif circuit.num_qubits == 5: 
    #     sentenceDM = partial_trace(state=rho, qargs=[0, 1, 2, 4])

    # # simple sentence like man prepares food, without any adjectives, man will be daggered there 
    # elif circuit.num_qubits == 3: 
    #     sentenceDM = partial_trace(state=rho, qargs=[0,2])



    '''
    I commented this out!!! 
    '''

    print('debug, shape density matrix: ', rho)

    print('Debug Debug', traceList)
    print(circuit)

    # if not addedQubit:
    #     sentenceDM = partial_trace(state=rho, qargs= traceList)
    # if addedQubit: 
    #     traceList = [i+1 for i in traceList]
    #     traceList.insert(0,0)
    #     # print('modified tracelist if appended:', traceList)
    #     sentenceDM = partial_trace(state=rho, qargs= traceList)


    '''
    I commented this out!! 
    '''


    # tracing out... 
    tracedOutRho = partial_trace(state = rho, qargs= qubitToTraceOut)

    print('we measure traceList before tracing out: ', traceList)

    # create new tracelist, because it has shifted due to the tracing out
    for i in range(len(traceList)): 
        if traceList[i] > qubitToTraceOut[0]:
            traceList[i] = traceList[i]-1

    print('we traced out ', qubitToTraceOut, ' and now we measure the traceList ', traceList)


    measureList = tracedOutRho.measure(traceList)


    print('the list of measurements or whatever: ', measureList)
    # while '1' in measureList[0]: 
    #     # print('how often does this repeat itself?')
    #     print(measureList[0])
    #     print('iteration')
    #     measureList = tracedOutRho.measure(traceList)


    sentenceDM = measureList[1]

    print('sentence DM is, after tracing out:')
    print(sentenceDM)

    # try tracing out all other qubits after measuring 
    # traceList.append(qubitToTraceOut[0])

    print('we trace out: ', traceList)

   


    # let's take a small detour to compare the resulting density matrix to the correct pure density matrix of the corresponding category 
    from qiskit.quantum_info import state_fidelity

    iTMatrix = np.array([[1,0], [0,0]])
    foodMatrix = np.array([[0,0], [0,1]])
    iTDM = DensityMatrix(iTMatrix)
    foodDM = DensityMatrix(foodMatrix)

    
    # the dm that captures the meaning of the sentence on only the sentence qubit (making it 2-dim)
    sentenceDMOneQubit = partial_trace(sentenceDM, traceList)

    listOfFidelities = [state_fidelity(sentenceDMOneQubit, foodDM) , state_fidelity(sentenceDMOneQubit, iTDM)]

    

    print('sentence dm one qubit is: ', sentenceDMOneQubit)

    # else: 
    #     print('number of qubits: ', circuit.num_qubits)
    #     print(circuit)

    return sentenceDMOneQubit, listOfFidelities

'''
now for the application of the code
'''

def Main(listOfCircuits: list, parameterDict: dict, wordsToForget: list, amplitudeEncoded: bool, traceOutInsteadOfMeasure: bool = False):
    '''
    Iterate over all circuits in the test_circuits list and for each of them create the density matrix as explained above 
    ''' 

    listOfDensityMatrices: list = []


    from pytket.extensions.qiskit import tk_to_qiskit

    numQubitsList: list = []

    for i in range(len(listOfCircuits)): 

        print('iteration number: ', i)

        tempCirc = listOfCircuits[i]

        # if i == 8: 
        #     tempCirc.draw()

        tempCircTK = tempCirc.to_tk()

        tempCircQiskit = tk_to_qiskit(tempCircTK)

        # if i == 8: 
        #     print(tempCircQiskit)

        # circuitData = tempCircQiskit.data

        

        alteredTempCirc, addedQubit, traceList, qubitToTraceOut = SetUpDiags(tempCircQiskit, parameterDict, amplitudeEncoded, wordsToForget, traceOutInsteadOfMeasure)

        


        # tracelist is list of qubits that we trace over. These are the qubits that are measured in the circuit 
        # list of qubits that shall be forgotten when traced over
        # traceList = []
        # circuitData = alteredTempCirc.data
        # for i in range(len(circuitData)):
        #     if str(circuitData[i].operation.name) == 'measure':
        #         traceList.append(circuitData[i].qubits[0].index)
        # print('traceList: ', traceList)


        # to access density matrix later on ONLY FOR AER SIMULATOR
        # alteredTempCirc.save_density_matrix(conditional = True) 


        densityMatrix, dmOneQubit = GetDensityMatrix(alteredTempCirc, traceList, addedQubit, traceOutInsteadOfMeasure, qubitToTraceOut) 
        
        listOfDensityMatrices.append(densityMatrix)

        numQubitsList.append(tempCircQiskit.num_qubits)

    return listOfDensityMatrices, numQubitsList, dmOneQubit
    



