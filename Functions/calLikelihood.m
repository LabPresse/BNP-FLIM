function LikeOut = calLikelihood(Data,EmpParam,S,Lambda,Ntmp,T_IRF,T,Sig_IRF)
%calLikelihood finds the likelihood of the proposed lamda
%
%INPUT:
%   Data: Structure array containing photon arrival times (ns)
%   EmpParam: Structure containing parameters of the experiment
%   S:  Sampled S in the previous iteration
%   Lambda: Proposed lambda (1/ns)
%
%OUTPUT:
%   Likelihood: Calculated likelihood of the proposed lambda
%
%Created by:
%   Mohamadreza Fazel (Presse Lab, 2020)
%

if length(Lambda)==1
    LambdaS = Lambda*ones(length(S),1);
else
    LambdaS = Lambda(S)';
end

LambdaT = repmat(LambdaS,[1,Ntmp+1]);
DataT = repmat(Data,[1,Ntmp+1]);
NT = repmat((0:Ntmp),[length(Data),1]);
  

LikeExp = (LambdaT/2).*exp((LambdaT/2).*(2*(T_IRF-DataT-NT*T) + ...
        LambdaT*Sig_IRF.^2));
LikeErf = erfc((T_IRF-DataT-NT*T+LambdaT*Sig_IRF.^2) ...
        /(sqrt(2)*Sig_IRF)); 
    
LikeOut = sum(LikeExp.*LikeErf,2);
LikeOut = sum(log(LikeOut));

end
