function[Lambda,LikeL] = sampleLambda_Unif(Data,EmpParam,S,Lambda_Old,Alpha,Beta,Alpha_Prop,Ntmp,T_IRF,T_intv,Sig_IRF)
%sampleLambda uses a MH algorithm to sample lambda
%
%INPUT:
%   Data: Structure array containing photon arrival times (ns)
%   EmpParam: Structure containing parameters of the experiment
%   S: Sampled S parameters in the previous iteration
%   Lambda_Old: Sampled lambda in the previous iteration
%   Alpha: Shape parameter of gamma prior on lambda
%   Beta:  Rate parameter of gamma prior on lambda
%   Alpha_Prop: Shape parameter for lambda proposal distribution
%
%OUTPUT:
%   Lambda: Sampled lambda
%
%Created by:
%   Mohamadreza Fazel (Presse Lab, 2020) 
%

Stmp = S.S(Data.W==1);
Dt = Data.Dt;

Like_Old = calLikelihood(Dt,EmpParam,Stmp,Lambda_Old,Ntmp,T_IRF,T_intv,Sig_IRF);
Lambda_Prop = gamrnd(Alpha_Prop,Lambda_Old/Alpha_Prop);

Like_Prop = calLikelihood(Dt,EmpParam,Stmp,Lambda_Prop,Ntmp,T_IRF,T_intv,Sig_IRF);
LikeRatio = sum(Like_Prop-Like_Old);
PriorRatio = sum(log(gampdf(Lambda_Prop,Alpha,Beta))-log(gampdf(Lambda_Old,Alpha,Beta)));
PropRatio = sum(log(gampdf(Lambda_Old,Alpha_Prop,Lambda_Prop/Alpha_Prop))...
     -log(gampdf(Lambda_Old,Alpha_Prop,Lambda_Old/Alpha_Prop)));
 
A = LikeRatio+PriorRatio+PropRatio;
if A > log(rand())
    Lambda = Lambda_Prop;
    LikeL = Like_Prop;
else
    Lambda = Lambda_Old;
    LikeL = Like_Old;
end

end
