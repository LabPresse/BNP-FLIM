function [Rho_IndP,Rho_Found,AcceptGP]=sampleRho(Data,EmpParam,Chol_IndP,...
    S,Rho_IndP,Rho_Found,D,TestP,IDTestP,AcceptGP,Matrix,Sig_GP,MeanXi,M,DEBUG)
%sampleRho jointly samples the Rho-profiles using MH method
%
%INPUTS: 
%   Data:      Input data
%   EmpParam:  Structure containing parameters required to simulate the data
%   Chol_IndP: Cholesky decomposition of the correlation between inducing inducing points
%   S:         Set of indicator parameters
%   Rho_IndP:  Learned profiles at inducing points
%   Rho_Found: Learned profiles at test points
%   D:         Size of the grid in grid of test points
%   TestP:     Coordinates of the test points grid
%   IDTestP:   Matrix of flags indicating what test points are no further
%              than 3sigma from each pixel.
%   AcceptGP:  Number of accepted proposed Rho profiles so far
%   Matrix:    Correlation matrix between test and inducing points
%   Sig_GP:    Mean size of proposed jumps in profiles
%   MeanXi:    Mean of the GP priors
%   M:         Number of lifetime species
%   DENUG:     If 1 shows an animation of the current profiles (This is slow)
%
%OUTPUT:
%   Rho_IndP:  Current profiles at inducing points
%   Rho_Found: Current profiles at test points
%   AcceptGP:  Updated accepted number of proposed profiles
%
%Created by:
%   Mohamadreza Fazel (Presse lab 2020)
%


%Updating first Rho
Rho_IndPtmp = Rho_IndP;
Sig = Sig_GP; %Sigma of drawn random number used in sampling from multivariate normal
Rho_tmp = Rho_Found;

for mm = 1:M
    
    Xi = log(Rho_IndP(mm).Rho);
    tXi = Xi + Sig*((randn([1,length(Rho_IndP(1).Rho)])*Chol_IndP)');
    Rho_IndPtmp(mm).Rho = exp(tXi);
    Rho_tmp(mm).Rho = exp(Matrix*tXi);        

end

NormFlag = 1;
LogLike_Old = 0;
for ii = 1:length(Data)
     [Pi_S,P0A] = calPi(EmpParam,Rho_Found,TestP,D,ii,IDTestP,NormFlag);
     St = S(ii).S;
     LogLike_Old = LogLike_Old + sum(log(Pi_S(St))) + ...
       sum(Data(ii).W==1)*log(1-prod(P0A)) + ...
       sum(Data(ii).W==0)*(sum(log(P0A)));
end
    
LogLike_Prop = 0;
for ii = 1:length(Data)
     [Pi_S,P0A] = calPi(EmpParam,Rho_tmp,TestP,D,ii,IDTestP,NormFlag);
     St = S(ii).S;
     LogLike_Prop = LogLike_Prop + sum(log(Pi_S(St))) + ...
         sum(Data(ii).W==1)*log(1-prod(P0A)) + ...
         sum(Data(ii).W==0)*(sum(log(P0A)));
end

LogPriorProp = 0;
LogPriorOld = 0;
for mm = 1:M
    LogPriorProp = LogPriorProp -0.5*sum(((log(Rho_IndPtmp(mm).Rho)-(MeanXi(mm)))'/Chol_IndP).^2) ...
        - sum(log(diag(Chol_IndP))) - length(Rho_IndPtmp(mm).Rho)*log(2*pi)/2;
    LogPriorOld = LogPriorOld -0.5*sum(((log(Rho_IndP(mm).Rho)-(MeanXi(mm)))'/Chol_IndP).^2) ...
        - sum(log(diag(Chol_IndP))) - length(Rho_IndP(mm).Rho)*log(2*pi)/2;
end
LogPriorRatio = LogPriorProp - LogPriorOld;

if -abs(LogLike_Prop) == LogLike_Prop
   if LogLike_Prop - LogLike_Old + LogPriorRatio > log(rand())
       Rho_Found = Rho_tmp; 
       Rho_IndP = Rho_IndPtmp;
       AcceptGP = AcceptGP + 1;
   end
end

for mm = 1:M
    Rho_Found(mm).Rho = single(Rho_Found(mm).Rho);
end

if DEBUG
    figure(111);
    L = [max(EmpParam.Xi_Y),max(EmpParam.Xi_X)]/EmpParam.PixelSize + 0.5;
    [Xg,Yg] = meshgrid((0.5:L(2)-0.5),(0.5:L(1)-0.5));
    for mm = 1:M
        Rho = Rho_IndP(mm).Rho;
        Rho1 = reshape(Rho,L);
        surf(Xg,Yg,Rho1)
        if mm == 1
            hold(gca,'on');
        end
    end
    xlabel('X');ylabel('Y')
    pause(0.02)
    hold off;
end

end
