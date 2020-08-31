function Data = tomography_pauli(MaxNumberOfQubits,NumberOfTrials)
    % runtime takes as inputs the maximum number of qubits MaxNumOfQubits 
    % and the maximum number of copies of the state (in each case
    % dimension) MaxCopies. The function returns "Data", a matrix whose
    % first row gives the number of qubits, the second row is the runtime
    % taken by qse_apg to reconstruct the state and the third row is the
    % infidelity between the Generator state and the resulting state
    % found by qse_apg
    
   
    Data = zeros(3, MaxNumberOfQubits);
    
    for NumOfQubits=1:MaxNumberOfQubits
        
        d = 2^NumOfQubits;
        
        %v = RandomStateVector(d);
        
        Proj = {};

            for n=1:NumOfQubits

                Proj = cat(2,Projectors(1), Proj);

            end
        
        copies = 3^NumOfQubits*500*d;
        [Average_runtime, Average_fidelity] = experimentSimulator(NumberOfTrials,d,copies,Proj);
        
        fprintf('%2d    %2.5f    %2.5f\n', [NumOfQubits , Average_runtime, Average_fidelity]')
        
        Data(:,NumOfQubits) = [NumOfQubits ; Average_runtime ; Average_fidelity ];
        
        
    end
    
    name = 'data_Paulis.mat';
        
    save(name ,'Data');
    
    
end
