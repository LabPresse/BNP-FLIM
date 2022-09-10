function [LogP,Like_Dt] = calPost_Unif(Data,Chain,S,EmpParam,BNP,Ntmp)

Like_Dt = 0;
D = BNP.D;
Lambda = Chain.Lambda;
for ii = 1:length(Data)
     Stmp = S(ii).S(Data(ii).W==1);
     Like_Dt = Like_Dt + calLikelihood(Data(ii).Dt,EmpParam,Stmp,Lambda,Ntmp,...
         EmpParam.T_IRF,EmpParam.T,EmpParam.Sig_IRF);
end

M = BNP.M;
P0A = zeros(1,M);
Pi_S = zeros(1,M);
Rho_IndP = Chain.Rho_IndP;
for mm = 1:M
    P0A(mm) = exp(-EmpParam.Mu(mm)*EmpParam.Dp*EmpParam.OmegaX^2*EmpParam.OmegaZ*pi*Rho_IndP(mm).Rho/2);
end
for mm = 1:M
    Pi_S(mm) = (1-P0A(mm))*prod(P0A)/P0A(mm);
end
Pi_S = Pi_S/sum(Pi_S);

Like_Rho = 0;
for ii = 1:length(Data)
    St = S(ii).S;
    Like_Rho = Like_Rho + sum(log(Pi_S(St))) + ...
        sum(Data(ii).W==1)*log(1-prod(P0A)) + ...
        sum(Data(ii).W==0)*sum(log(P0A));
        
end
    
PriorGamma = sum(log(gampdf(Lambda,BNP.Alpha,BNP.Beta)));
PriorRho = 0;
for ii = 1:length(Chain.Rho)
    tRho = Chain.Rho_IndP(ii).Rho;
    PriorRho = PriorRho + log(gampdf(tRho,1,20));
end
LogP = Like_Dt + Like_Rho + PriorGamma + PriorRho;

end
