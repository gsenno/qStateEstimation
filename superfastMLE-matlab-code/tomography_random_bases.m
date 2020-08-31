function Data = tomography_random_bases( MaxNumOfQubits, NoT, mubs)
    % runtime takes as inputs the maximum number of qubits MaxNumOfQubits 
    % and the maximum number of copies of the state (in each case
    % dimension) MaxCopies. The function returns "Data", a matrix whose
    % first row gives the number of qubits, the second row is the runtime
    % taken by qse_apg to reconstruct the state and the third row is the
    % infidelity between the Generator state and the resulting state
    % found by qse_apg
    
    NoC = 60000;
    
    Data = zeros(3, MaxNumOfQubits);
    
    
    for NumOfQubits=1:MaxNumOfQubits
        
        d = 2^NumOfQubits;    
        ruido = 0.1;
        
        bases=zeros(d,d*(d+1));
        for i=0:d
            bases(:,i*d+1:(i+1)*d) = RandomUnitary(d);
        end
        cumulative_time = 0;
        cumulative_fidelity = 0;
        
        % NoT: Number of trials
        for trials=1:NoT
            
            %rho_target = (v*v');

            rho_target = RandomDensityMatrix(d);

            noisy_rho = (1-ruido)*rho_target + ruido*eye(d)/d;
            
            probabilidades0 =  qmt( noisy_rho ,  bases);
            norma = sum(probabilidades0);
            probabilidades1 = probabilidades0/norma;

            copies = (d+1)*100*d;
            counts = histc(rand(copies,1), [0; cumsum(probabilidades1)]);
            probabilidades2 = norma*counts(1:end-1)/(1000*d);
                
             
            tic
            
            result = qse_apg(bases, probabilidades2);  % Reconstruction of the state using Proj
            
            cumulative_time = cumulative_time + toc;
            
            %QETLABs' Fidelity is the sqroot of the more standard fidelity
            fidelity = Fidelity( rho_target, result )^2;
            
            cumulative_fidelity = cumulative_fidelity + fidelity;

        end
        
        
        Average_runtime = cumulative_time/NoT;
        
        Average_fidelity = cumulative_fidelity/NoT;
        
        fprintf('%2d    %2.5f    %2.5f\n', [NumOfQubits , Average_runtime, Average_fidelity]')
        
        Data(:,NumOfQubits) = [NumOfQubits ; Average_runtime ; Average_fidelity ];
        
        
    end 
    
    
    %Date = datetime(now,'ConvertFrom','datenum');
       
    dir = pwd;
        
    name = 'data_rand.mat';
        
    %save(dir + Date + name ,'Data');
    save(name ,'Data');
    
end

% temp = struct2cell( load('MUBs1.mat') );
% mubs1 = horzcat(temp{:});    
% temp = struct2cell( load('MUBs2.mat') );
% mubs2 = horzcat(temp{:});  
% temp = struct2cell( load('MUBs3.mat') );
% mubs3 = horzcat(temp{:});  
% temp = struct2cell( load('MUBs4.mat') );
% mubs4 = horzcat(temp{:});  
% temp = struct2cell( load('MUBs5.mat') );
% mubs5 = horzcat(temp{:});  
% temp = struct2cell( load('MUBs6.mat') );
% mubs6 = horzcat(temp{:});  
% temp = struct2cell( load('MUBs7.mat') );
% mubs7 = horzcat(temp{:}); 
% temp = struct2cell( load('MUBs8.mat') );
% mubs8 = horzcat(temp{:});
% 
%  mubs{1} = mubs1;
%  mubs{2} = mubs2;
%  mubs{3} = mubs3;
%  mubs{4} = mubs4;
%  mubs{5} = mubs5;
%  mubs{6} = mubs6;
%  mubs{7} = mubs7;
%  mubs{8} = mubs8;
% 
