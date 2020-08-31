import numpy as np
import qutip as qt

import itertools as it
import matplotlib.pyplot as plt
import functools as ft
from quantum_state_estimation import TomographyExperiment, qStateEstimation

def experimentSimulator(dimension,measurementType,measurements,no_trials,copiesPerMeasurement):
    time_pio , fidelities_pio = 0 , 0
        
    for _ in range(no_trials):
        
        rho_target = qt.rand_dm(dimension).full()
        white_noise_level = 0.1
        rho = (1 - white_noise_level)*rho_target + white_noise_level*np.identity(dimension)/dimension
        
        probabilities = probability_distributions(measurements, dimension, rho, copiesPerMeasurement)
        experiment = TomographyExperiment(measurementType,zip(measurements,probabilities))
        
        rho_estimated , runtime = qStateEstimation(experiment, 
                                        dimension = dimension , 
                                        min_HS_distance = 10e-6,
                                        max_iterations = 1
                                        )
        
        fidelity = qt.fidelity( qt.Qobj(rho_estimated) , qt.Qobj(rho_target) )
        
        time_pio = time_pio + runtime
        fidelities_pio = fidelities_pio + fidelity
    
    return time_pio/no_trials,fidelities_pio/no_trials

def born_rule_probabilities(bases, rho, numberOfBases,dimension):
    probabilities=np.zeros((numberOfBases,dimension))
    for i,basis in enumerate(bases):
        probabilities[i]=np.real( np.diag(np.asmatrix( basis ).getH() @ rho @ basis ))
    return probabilities

def probability_distributions(bases, dimension, rho, no_countings ):
    
    probabilities = np.zeros((len(bases),dimension))
    for i,basis in enumerate(bases):     
         
        pdi = np.real( np.diag(np.asmatrix( basis ).getH() @ rho @ basis ) )
        pdi  = np.array( [0] + list(pdi) )
        pdi, _ = np.histogram(list(np.random.random_sample(no_countings)), np.cumsum(pdi))

        probabilities[i] = pdi/np.sum(pdi)
    
    return probabilities

def pauli_bases(num_of_qubits):
    
    basisZ = np.matrix([[1,0],[0,1]])
    basisX = np.matrix([[1,1],[1,-1]])/np.sqrt(2)
    basisY = np.matrix([[1,1],[1j,-1j]])/np.sqrt(2)
    
    qubitBasis = [basisZ,basisX,basisY]
    
    if num_of_qubits == 1:
        return qubitBasis
    else:
        pauliIndexes = list(it.product(range(3),repeat=num_of_qubits))
        return [ft.reduce(lambda x, y: np.kron(x,y),[qubitBasis[i] for i in index]) for index in pauliIndexes]


class PauliBases:
    basisZ = np.matrix([[1,0],[0,1]])
    basisX = np.matrix([[1,1],[1,-1]])/np.sqrt(2)
    basisY = np.matrix([[1,1],[1j,-1j]])/np.sqrt(2)
    
    qubitBases = [basisZ,basisX,basisY]
    
    no_qubits = 2
    
    def __init__(self,numberOfQubits):
        self.no_qubits=numberOfQubits
    
    def _generatorRecursion(self,num_of_qubits):
        if num_of_qubits == 1:
            yield from self.qubitBases
        else:
            yield from (np.einsum('ik,jl', prod, aBasis).reshape(2**num_of_qubits,2**num_of_qubits)
                         for aBasis in self.qubitBases 
                         for prod in self._generatorRecursion(num_of_qubits-1))
    
    def generator(self):
        return self._generatorRecursion(self.no_qubits)
        
import matplotlib.pyplot as plt

def plot_data( time_pio , fidelity_pio, runtime_mle, fidelity_mle, name_of_file = 'data.png' ):
    
    
    min_noq = list( time_pio.keys() )[0]
    
    max_noq = list( time_pio.keys() )[-1]

    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(17,6))

    plt.rc('xtick', labelsize= 14) 

    plt.rc('ytick', labelsize= 14) 

    #fig.tight_layout()

    #plt.xlabel("Number of Qubits", fontsize=14)

    ######### runtime ########################################################
    
    ax1.plot( list( runtime_mle.keys() ), list( runtime_mle.values() ) , linewidth = 2.4, markersize = 9 , 
             linestyle='-', mfc='none', marker='^', color='g')

    ax1.plot( list( time_pio.keys() ), list( time_pio.values() ), linewidth = 2.4, markersize = 9 , 
             linestyle='--', mfc='none', marker='o', color='m')

    ax1.set_yscale('log')

    ax1.set_ylabel("t (s)", fontsize = 15)

    ax1.set_xlabel("Number of Qubits", fontsize = 15)

    ax1.set_xticks( list( time_pio.keys() ) )

    #ax1.set_title("Runtime", fontsize = 16)

    ax1.yaxis.set_ticks_position('both')

    ax1.xaxis.set_ticks_position('both')

    ax1.tick_params(which='major', direction='in', length=6, width=1.5)

    ax1.tick_params(which='minor', direction='in', length=2.5, width=1.5)

    ax1.legend( ['Superfast MLE', 'Algorithm 1'], fontsize = 14, loc = 'lower right')

    ax1.set_xlim( min_noq , max_noq )
    
    #ax.tick_params(which='minor', length=4, color='r')
    
    #ax1.tick_params(labeltop=True, labelright=True)


    ######### fidelity ########################################################
    
#     ax2.plot( list( fidelity_pio.keys() ), list( fidelity_pio.values() ),
#              linewidth = 2.4, markersize = 6,  linestyle='-', marker='s', color='b')
    
#     ax2.plot( list( fidelity_pio.keys() ), list( fidelity_pio.values() ),
#              linewidth = 2.4, markersize = 6,  linestyle='-', marker='s', color='b')


    ax2.plot( list( fidelity_mle.keys() ), list( fidelity_mle.values() ) , linewidth = 2.4, markersize = 9 , 
             linestyle='-',  marker='^', color='g')

    ax2.plot( list( fidelity_pio.keys() ), list( fidelity_pio.values() ),
             linewidth = 2.4, markersize = 9 , linestyle='--',  marker='o', color='m')

    #ax2.set_title('Fidelity', fontsize=16)

    ax2.set_ylabel("Fidelity", fontsize = 15)

    ax2.set_xlabel("Number of Qubits", fontsize = 15)

    ax2.set_xticks( list( fidelity_pio.keys() ) )

    ax2.yaxis.set_ticks_position('both')

    ax2.xaxis.set_ticks_position('both')

    ax2.tick_params( direction='in', length = 6, width = 1.5 )

    ax2.set_ylim(0,1.01)
    
    ax2.legend( ['Superfast MLE', 'Algorithm 1'], fontsize = 14, loc = 'lower right')

    ax2.set_xlim( min_noq , max_noq )

    plt.tight_layout()

    plt.savefig( name_of_file )
    

    return plt.show()

