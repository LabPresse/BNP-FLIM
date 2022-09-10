function [Rho_IndP,AcceptGP]=sampleRho_Unif(Data,EmpParam,S,Rho_IndP,AcceptGP)
%Updating first Rho
Rho_IndPtmp = Rho_IndP;
Sig = 25000; %Sigma of drawn random number used in sampling from multivariate normal
M = length(Rho_IndP);

for mm = 1:M
    Rho_IndPtmp(mm).Rho = gamrnd(Sig,Rho_IndP(mm).Rho/Sig);    
end
       
P0A = zeros(1,M);
Pi_S = zeros(1,M);
for mm = 1:M
    P0A(mm) = exp(-EmpParam.Mu(mm)*EmpParam.Dp*EmpParam.OmegaX^2*EmpParam.OmegaZ*pi*Rho_IndP(mm).Rho/2);
end
for mm = 1:M
    Pi_S(mm) = (1-P0A(mm))*prod(P0A)/P0A(mm);
end
Pi_S = Pi_S/sum(Pi_S);

LogLike_Old = 0;
for ii = 1:length(Data)
    St = S(ii).S;
    LogLike_Old = LogLike_Old + sum(log(Pi_S(St))) + ...
        sum(Data(ii).W==1)*log(1-prod(P0A)) + ...
        sum(Data(ii).W==0)*sum(log(P0A));
        
end
   
P0A = zeros(1,M);
Pi_S = zeros(1,M);
for mm = 1:M
    P0A(mm) = exp(-EmpParam.Mu(mm)*EmpParam.Dp*EmpParam.OmegaX^2*EmpParam.OmegaZ*pi*Rho_IndPtmp(mm).Rho/2);
end
for mm = 1:M
    Pi_S(mm) = (1-P0A(mm))*prod(P0A)/P0A(mm);
end
Pi_S = Pi_S/sum(Pi_S);
LogLike_Prop = 0;
for ii = 1:length(Data)
     St = S(ii).S;
     LogLike_Prop = LogLike_Prop + sum(log(Pi_S(St))) + ...
         sum(Data(ii).W==1)*log(1-prod(P0A)) + ...
         sum(Data(ii).W==0)*sum(log(P0A));
end

LogLikeR = LogLike_Prop - LogLike_Old;
LogPriorR = 0;
LogPropR = 0;
for mm = 1:M
    LogPriorR = LogPriorR+log(gampdf(Rho_IndPtmp(mm).Rho,1.01,150))-log(gampdf(Rho_IndP(mm).Rho,1.01,150));
    LogPropR = LogPropR+log(gampdf(Rho_IndP(mm).Rho,Sig,Rho_IndPtmp(mm).Rho/Sig))...
        -log(gampdf(Rho_IndPtmp(mm).Rho,Sig,Rho_IndP(mm).Rho/Sig));
end

if  LogLikeR+LogPriorR+LogPropR > log(rand())
    Rho_IndP = Rho_IndPtmp;
    AcceptGP = AcceptGP + 1;
end

end
