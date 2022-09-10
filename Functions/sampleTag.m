function S = sampleTag(Data,EmpParam,Lambda,Rho_Found,TestP,IDTestP,D,Ntmp)
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

S(length(Data)).S = 0;
for ii = 1:length(Data)

    NormFlag = 1;
    [Pi_S,~] = calPi(EmpParam,Rho_Found,TestP,D,ii,IDTestP,NormFlag);

    L = length(Data(ii).W);
    Rnd1 = rand(L,1);
    Pi_S = Pi_S';
    [~,St] = max(cumsum(Pi_S,2)-Rnd1>0,[],2);
    S(ii).S = St;
    
    Dt = Data(ii).Dt;
    [Probs,~] = calProbs(Dt,EmpParam,Lambda,Pi_S,Ntmp);
    Rnd2 = rand(length(Data(ii).Dt),1);
    [~,Stmp] = max(cumsum(Probs,2)-Rnd2>0,[],2);
    
    S(ii).S(Data(ii).W==1) = Stmp;
        
end

end
