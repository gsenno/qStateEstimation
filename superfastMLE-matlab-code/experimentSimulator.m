function [Average_runtime,Average_fidelity] = experimentSimulator(NumberOfTrials,dimension,copies,measurements)

    cumulative_time = 0;
        
    cumulative_fidelity = 0;
        
    for trial=1:NumberOfTrials

        ruido = 0.1;

        %rho_target = (v*v');

        rho_target = RandomDensityMatrix(dimension);

        noisy_rho = (1-ruido)*rho_target + ruido*eye(dimension)/dimension;



        probabilidades0 =  qmt( noisy_rho, measurements);

        norma = sum(probabilidades0);

        probabilidades1 = probabilidades0/norma;

        counts = histcounts(rand(copies,1), [0; cumsum(probabilidades1)]);
        counts = counts';
        probabilidades2 = norma*counts/copies;
        
        tic
        result = qse_apg(measurements, probabilidades2);  % Reconstruction of the state using Proj

        cumulative_time = cumulative_time + toc;

        fidelity = Fidelity( rho_target, result )^2;

        cumulative_fidelity = cumulative_fidelity + fidelity;

    end


    Average_runtime = cumulative_time/NumberOfTrials;

    Average_fidelity = cumulative_fidelity/NumberOfTrials;