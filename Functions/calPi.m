function [P,P0A] = calPi(EmpParam,RhoIn,TestP,D,ii,IDTestP,NormalizeFlag)
%calPi finds the probability of photons from different species given concentrations
%
%INPUT:
%   EmpParam: Structure containing parameters of the experiment
%   RhoIn: Structure containing inline functions of concentrations
%   ii: Index of the confocal spot
%
%OUTPUT:
%   P: Vector of calculated probabilities
%
%Created by:
%   Mohamadreza Fazel (Presse Lab, 2020)
%   

M = length(RhoIn);
P0A = zeros(M,1);
Pi_At = zeros(M,1);
for mm = 1:M
    IntA = discreteInt(TestP,RhoIn(mm).Rho,EmpParam,D,ii,IDTestP);
    P0A(mm) = exp(-EmpParam.Mu(mm)*EmpParam.Dp*IntA);
end

for mm = 1:M
    Pi_At(mm) = (1-P0A(mm))*prod(P0A)/P0A(mm);
end

P=zeros(M,1);
if NormalizeFlag
    for mm = 1:M
        P(mm) = Pi_At(mm)/sum(Pi_At);
    end
else
    for mm = 1:M
        P(mm) = Pi_At(mm);
    end
end

end

