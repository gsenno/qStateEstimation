import numpy as np
import time
from scipy.optimize import root
import enum

class TomographyExperiment:
    def __init__(self,measurementType,measurementsAndProbabilities):
        self.measurementType = measurementType
        self.measurementsAndProbabilities = measurementsAndProbabilities
    
    def getMeasurementsAndProbabilities(self):
        return self.measurementsAndProbabilities
    
    def getMeasurementType(self):
        return self.measurementType

class MeasurementType(enum.Enum):
    Bases = 0
    PVMs = 1
    GeneralPOVMs = 2
    
def T(measurementType,rho_0, measurement , probabilities ):
    if measurementType==MeasurementType.Bases:
        aBasis=measurement
        rotatedState = np.asmatrix( aBasis ).getH() @ rho_0 @ aBasis 
        stateWithImposedInformation = rotatedState - np.diag( np.diag( rotatedState ) ) + np.real( np.diag( probabilities ) )
        rho_1 = aBasis @ stateWithImposedInformation @ np.asmatrix( aBasis ).getH()
    else:
        if measurementType==MeasurementType.PVMs:
            for i,projector in enumerate(measurement):
                rho_1 = rho_0 + (
                    probabilities[i]-np.trace(projector@rho_0)
                    )/np.trace(projector).real*projector
            else:
                rho_aux = rho_0
                for i,effect in enumerate(measurement):
                    rho_1 = rho_aux + (probabilities[i]-np.trace(effect@rho_aux)
                                     )/np.trace(effect@effect).real*effect
                    rho_aux = rho_1      
    return rho_1

def applyCompositionT(measurementType,rho_0,measurementsAndProbabilities):
    rho = rho_0
    for measurement,probabilities in measurementsAndProbabilities:
        rho = T(measurementType,rho_0, measurement, probabilities )
        rho_0 = rho
    return rho

def projectToDensityOperator(rho,dimension):
    
    eigenvalues , eigenvectors = np.linalg.eigh( rho )
    eigenvalues = np.real(eigenvalues)
    U = np.matrix(eigenvectors)
    
    f = lambda x : np.sum([abs(eigval-x) for eigval in eigenvalues]) - dimension*x - np.trace(rho).real
    x0=root(f,0).x[0].real
    
    projection = U@np.diag([max([eigval-x0,0]) for eigval in eigenvalues])@U.getH()
    
    return projection/(np.trace(projection).real)

def qStateEstimation(tomographyExperiment, dimension = 2, min_HS_distance = 1e-10, max_iterations = 1e10):
    
    start_time_pio = time.time()
    
    #initial estimate
    seed = np.eye(dimension)
    rho = applyCompositionT(tomographyExperiment.getMeasurementType(),
                            seed,
                            tomographyExperiment.getMeasurementsAndProbabilities())
    
    #iteration until convergence in HS distance
    hs_distance = 2.0
    n=1
    while (hs_distance > min_HS_distance) & (n < max_iterations):
        rho_0 = rho
        rho = applyCompositionT(tomographyExperiment.getMeasurementType(),
                                rho_0,
                                tomographyExperiment.getMeasurementsAndProbabilities())
        hs_distance = HS_distance(rho,rho_0)
        n += 1
        
    #find the density operator closest to rho in Frobenius norm.
    rho_estimated = projectToDensityOperator(rho,dimension)
    
    dt_pio = time.time() - start_time_pio

    return rho_estimated, dt_pio


def HS_distance(rho,sigma):
    return np.trace(np.linalg.matrix_power((rho - sigma),2))
