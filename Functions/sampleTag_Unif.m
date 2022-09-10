function S = sampleTag_Unif(Data,EmpParam,Lambda,Rho_IndP,Ntmp)
%sampleTag uses Gibbs algorithm to sample the parameter S
%
%INPUT:
%   Data: Structure array containing photon arrival times (ns)
%   EmpParam: Structure containing parameters of the experiment
%   Lambda: Sampled lambda in the previous iteration (1/ns)
%   Rho_Found: 
%   Pi_T: Relative probability of photons coming from different species
%   
%OUTPUT:
%   S: The sampled S parameters
%
%Created by:
%   Mohamadreza Fazel (Presse Lab, 2020)
%

M = length(Lambda);
P0A = zeros(1,M);
Pi_S = zeros(1,M);
for mm = 1:M
    P0A(mm) = exp(-EmpParam.Mu(mm)*EmpParam.Dp*EmpParam.OmegaX^2*EmpParam.OmegaZ*pi*Rho_IndP(mm).Rho/2);
end
for mm = 1:M
    Pi_S(mm) = (1-P0A(mm))*prod(P0A)/P0A(mm);
end
Pi_S = Pi_S/sum(Pi_S);

S(length(Data)).S = 0;
for ii = 1:length(Data)
    
    L = length(Data(ii).W);
    Rnd1 = rand(L,1);
    [~,St] = max(cumsum(Pi_S,2)-Rnd1>0,[],2);
    S(ii).S = St;
    
    Dt = Data(ii).Dt;
    [Probs,~] = calProbs(Dt,EmpParam,Lambda,Pi_S,Ntmp);
    Rnd2 = rand(length(Data(ii).Dt),1);
    [~,Stmp] = max(cumsum(Probs,2)-Rnd2>0,[],2);
    S(ii).S(Data(ii).W==1) = Stmp;
    
end

end



