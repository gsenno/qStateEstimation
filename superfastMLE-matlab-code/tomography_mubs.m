function Data = tomography_mubs( MaxNumOfQubits, NumberOfTrials)
    % runtime takes as inputs the maximum number of qubits MaxNumOfQubits 
    % and the maximum number of copies of the state (in each case
    % dimension) MaxCopies. The function returns "Data", a matrix whose
    % first row gives the number of qubits, the second row is the runtime
    % taken by qse_apg to reconstruct the state and the third row is the
    % infidelity between the Generator state and the resulting state
    % found by qse_apg
    
    Data = zeros(3, MaxNumOfQubits);
    
    temp = struct2cell( load('MUBs1.mat') );
    mubs{1} = horzcat(temp{:});    
    temp = struct2cell( load('MUBs2.mat') );
    mubs{2} = horzcat(temp{:});  
    temp = struct2cell( load('MUBs3.mat') );
    mubs{3} = horzcat(temp{:});  
    temp = struct2cell( load('MUBs4.mat') );
    mubs{4} = horzcat(temp{:});  
    temp = struct2cell( load('MUBs5.mat') );
    mubs{5} = horzcat(temp{:});  
    temp = struct2cell( load('MUBs6.mat') );
    mubs{6} = horzcat(temp{:});  
    temp = struct2cell( load('MUBs7.mat') );
    mubs{7} = horzcat(temp{:}); 
    temp = struct2cell( load('MUBs8.mat') );
    mubs{8} = horzcat(temp{:});
    
    for NumOfQubits=1:MaxNumOfQubits
        
        d = 2^NumOfQubits;    
        
        copies = (d+1)*100*d;
        [Average_runtime, Average_fidelity] = experimentSimulator(NumberOfTrials,d,copies,mubs{NumOfQubits});
        
        fprintf('%2d    %2.5f    %2.5f\n', [NumOfQubits , Average_runtime, Average_fidelity]')
        
        Data(:,NumOfQubits) = [NumOfQubits ; Average_runtime ; Average_fidelity ];
        
    end 
        
    name = 'data_mubs.mat';
        
    save(name ,'Data');
    
end
