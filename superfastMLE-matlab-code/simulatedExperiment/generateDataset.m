%This function generates a random density matrix 'rho' with a specified
%dimensionality and purity. It generates a measurement matrix with many
%options (see commented text below).
%input
%param: parameters of the gererated density matrix and measurement matrix.
%outputs:
%rho: randomly generated density matrix
%A: measurement matrix (see options below)
%data: data generated from rho and A with multinomial or poisson noise 


function [rho_target, data] = generateDataset(param,A,ruido)

d = param.d;
purity = param.purity;
counts = param.counts;

x = 1:d; %eigenvalue index
lambda=0;
purityTemp=0;
%Generate exponialy decreasing eigenvalues with the specified purity
while purityTemp<purity
    lambda = lambda + 0.001; %increase std until reach correct purity
    lam = exp(-lambda*x); %exponential distribution of eigenvalues
    lamb = lam/sum(lam);
    purityTemp = sum(lamb.^2);
end
lambi(1:d) = lamb; %eigenvalues

%Generate completely random density matrix with predefined eigenvalues
%rho_target = makeRandomDensityMatrix(lambi);
rho_target = RandomDensityMatrix(d);
%vec=RandomStateVector(d);
%rho_target=vec*vec';
rho = (1-ruido)*rho_target + ruido*eye(d)/d;

%calculate exact probabilities
temp0 = (conj(A).*((A)*rho));
dataExact = real(sum(temp0,2))*counts;

% if the total number of detector clicks is small enough, add multinomial
% noise, else add Gaussian-approximated Poisson noise. 
%Multinomial noise is too computationally intensive beyond a million clicks
%if ceil(counts/d*length(dataExact)) < 2E5 && license('test', 'Statistics_Toolbox')
%    data = mnrnd(ceil(counts/d*length(dataExact)),dataExact/sum(dataExact))';
%else
    data = dataExact + randn(length(dataExact),1).*sqrt(dataExact);
    data = ceil(data);
    data(data<0)=0;
%end

%calculate log likelihood of Gaussian-approximated Poisson noise. This
%should be equal to unity on average. It's very much like a chi^2 figure of
%merit. Uncomment to see the value
% log_L_GP = sum( ((dataExact-data)./sqrt(dataExact+eps))  .^2)/length(data)



