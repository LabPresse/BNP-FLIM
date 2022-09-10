function [Chain,S]=runBNPs_FLIM(Data,EmpParam,BNP,Rho_init,Lambda_init,S_init,Str)
%runBNPs_FLIM implements complete BNPs analysis of FLIM data.
%
%INPUT:
%   Data: Structure array containing photons arrival times where each element 
%         of the array holds arrival times for a corresponding pixel (ns)
%   EmpParams: Structure containing experimental parameters
%       PixelSize: pixel size (mu m)
%       Xi_X:      Confocal X center (mu m)
%       Yi_X:      Confocal Y center (mu m)
%       Zi_X:      Confocal Z center (mu m)
%       Mu:        Array of excitation rates of differrent species (1/ns)
%       Tau:       Array of lifetimes of the species (ns)
%       Dp:        Average width of laser pulses (ns)
%       Sig_IRF:   Sigma of IRF (ns)
%       T_IRF:     Mean of IRF (ns)
%       T:         Time interval between two consecurive laser pulses (ns)
%       OmegaX:    HFWM of confocal PSF along X and Y axis (mu m)
%       OmegaZ:    HFWM of confocal PSF along Z axis (mu m)   
%       PSF:       The combination of detection and illumination PSFs (inline function)   
%   BNP: Structure containing parameters used in the analysis
%       M:            Number of lifetime components (species)
%       Alpha:        Shape parameter of gamma prior on inverse of lifetimes
%       Beta:         Scale parameter of gamma prior on inverse of lifetimes
%       Alpha_Prop:   Parameter of proposal dist. of inverse of lifetimes
%       NJump:        Number of samples from the posterior
%       D:            Grid size of the test points (resolution of profiles) (mu m)
%       T:            Parameter of GP prior
%       L:            Parameter of GP prior (correlation) (mu m)
%       N:            Number of preceeding laser pulses to be considered
%       Sig_GP:       Average jump size in proposing new profiles
%       Sig_Xi:       Average jump size in proposing new Xi
%       Sig_Prior_Xi: Sigma of Gaussian prior on Xi 
%       DEBUG:        Show the current profiles (1)
%       UnifRho:      Flag if the concentrations are uniform or not (Default is 0)
%   Rho_init: Initial values for the profiles at the inducing points/pixel
%   centers (optional)
%   Lambda_init: Initial values for labmda (optional)
%   S_init: Initial values for tag parameter S (optional)
%   Str:    The name used to save the chain (Default = 'Chain')
%
%OUTPUT:
%   Chain: Structure array containing sampled parameters
%       Lambda: Chain of sampled inverse of lifetimes (Lambda)
%       Rho_IndP: Chain of sampled profiles at inducing points/pixel centers
%       Rho: Chain of sampled profiles at test points 
%       Xi: Chain of sampled means of GP priors
%
%Created by:
%   Mohamadreza Fazel (Presse Lab, 2021)
%

if nargin < 7
   Str = 'Chain'; 
end

NJump = BNP.NJump;
Alpha = BNP.Alpha;
Beta = BNP.Beta;
Alpha_Prop = BNP.Alpha_Prop;
D = BNP.D;
T = BNP.T;
L = BNP.L;
N = BNP.N;
M = BNP.M;
%Extension of the profiles in the edges (Default 2 pixels)
if isfield(BNP,'Exten')
    Exten = BNP.Exten;
else
    Exten = 2;
end
if isfield(BNP,'DEBUG')
    DEBUG = BNP.DEBUG;
else
    DEBUG = 0;
end
if isfield(BNP,'PerSample')
    PerSample = BNP.PerSample;
else
    PerSample = 50;
end
if isfield(BNP,'UnifRho')
    UnifRho = BNP.UnifRho;
else
    UnifRho = 0;
end

Chain(floor(NJump/PerSample)).Lambda = zeros(1,M);
Cent_Confocal = zeros(length(Data),3);
for ii = 1:length(Data)
    Cent_Confocal(ii,1) = Data(ii).X_Confocal;
    Cent_Confocal(ii,2) = Data(ii).Y_Confocal;
    Cent_Confocal(ii,3) = Data(ii).Z_Confocal;
end

%Initializing the chain based on the given inputs
if nargin >= 6
    if size(Lambda_init,1)>1
        Chain(1).Lambda = Lambda_init';
    else
        Chain(1).Lambda = Lambda_init;
    end
    Chain(1).Rho_IndP = Rho_init;
elseif nargin == 5
    S_init(length(Data)).S=0;
    for ii = 1:length(Data)
        S_init(ii).S = randi(M,length(Data(ii).W),1);
    end
    if size(Lambda_init,1)>1
        Chain(1).Lambda = Lambda_init';
    else
        Chain(1).Lambda = Lambda_init;
    end
    Chain(1).Rho_IndP = Rho_init;
elseif nargin == 4
    S_init(length(Data)).S=0;
    for ii = 1:length(Data)
        S_init(ii).S = randi(M,length(Data(ii).W),1);
    end
    Chain(1).Lambda = gamrnd(BNP.Alpha,BNP.Beta,[1 M]);
    Chain(1).Rho_IndP = Rho_init;
elseif nargin == 3
    S_init(length(Data)).S=0;
    for ii = 1:length(Data)
        S_init(ii).S = randi(M,length(Data(ii).W),1);
    end
    Chain(1).Lambda = gamrnd(BNP.Alpha,BNP.Beta,[1 M]);
    for mm = 1:M
        Chain(1).Rho_IndP(mm).Rho = 0.5+10*rand()*ones(size(Data));
    end
end

%Covariance matrix between inducing points
K_IndP = calCOV(Cent_Confocal,Cent_Confocal,T,L);
MinRange = -(Exten*EmpParam.PixelSize-D/2)*[1, 1]; 
MaxRange(1) = (Data(end).X_Confocal+Exten*EmpParam.PixelSize);
MaxRange(2) = (Data(end).Y_Confocal+Exten*EmpParam.PixelSize);
[Xg,Yg,Zg] = meshgrid(MinRange(1):D:MaxRange(1),MinRange(2):D:MaxRange(2),0);

TestP = [Xg(:),Yg(:),Zg(:)];
IDTestP = zeros([size(TestP,1),length(Data)]);
for ii = 1:size(Cent_Confocal,1)
     Dis = pdist2(TestP(:,1:2),Cent_Confocal(ii,1:2));
     IDTestP(Dis<3*EmpParam.OmegaX,ii) = ii;
end

if UnifRho == 0
    %Covariance matrix between inducing and test points
    K_test_IndP = calCOV(TestP,Cent_Confocal,T,L);
    %Choleski decomposition of covariance matrix
    Chol_IndP = cholcov(K_IndP+1000*eps*eye(size(K_IndP)));
    %The intermediate matrix that relate indusing points to test points
    Matrix = K_test_IndP/K_IndP;
    %Initializing the profiles at the test points based on the given values
    %for test points and the correlation matrix.
    for mm = 1:M
        Chain(1).Rho(mm).Rho = exp(Matrix*log(Chain(1).Rho_IndP(mm).Rho)); 
    end
else
    if length(Chain(1).Rho_IndP(1).Rho) > 1
        for mm = 1:M
            Chain(1).Rho_IndP(mm).Rho = Chain(1).Rho_IndP(mm).Rho(1);
        end
    end
    for mm = 1:M 
        Chain(1).Rho(mm).Rho = Chain(1).Rho_IndP(mm).Rho*ones(size(Cent_Confocal,1),1);
    end
end

AcceptGP = 0;
S = S_init;
T_IRF = EmpParam.T_IRF;
T_intv = EmpParam.T;
Sig_IRF = EmpParam.Sig_IRF;

if UnifRho
   Dtmp.Dt=[];
   Dtmp.N = 0;
   for ii = 1:numel(Data)
       Dtmp.Dt = cat(1,Dtmp.Dt,Data(ii).Dt);
       Dtmp.N = Dtmp.N + length(Data(ii).W);
   end
   S = [];
   Dtmp.W = zeros(Dtmp.N,1);
   Dtmp.W(1:length(Dtmp.Dt))=1;
   S(1).S = randi(M,[length(Dtmp.W),1]);
end
Chain(1).Xi = 0.5+10*rand([1,M]);
AcceptXi = 0;
Ind = 1;
ThisRho_IndP = Chain(1).Rho_IndP;
ThisRho = Chain(1).Rho;
ThisLambda = Chain(1).Lambda;
ThisXi = Chain(1).Xi;
for jj = 2:NJump
    
    if jj/5000 == floor(jj/5000)
        fprintf('Jump %d out of %d\n',jj,NJump);
        fprintf('Accepted GP jumps: %d\n',AcceptGP);
    end
    %Non-uniform profiles
    if UnifRho == 0
        
        %Sampling profiles
        [ThisRho_IndP,ThisRho,AcceptGP]=sampleRho(Data,EmpParam,Chol_IndP,...
            S,ThisRho_IndP,ThisRho,D,TestP,IDTestP,AcceptGP,Matrix,...
            BNP.Sig_GP,ThisXi,BNP.M,DEBUG);
        %Sampling mean of GPs
        [ThisXi,AcceptXi] = sampleXi(ThisRho_IndP,ThisXi,BNP,Chol_IndP,AcceptXi);
        %Sampling assignment parameter
        S = sampleTag(Data,EmpParam,ThisLambda,ThisRho,TestP,IDTestP,D,N);
        %Sampling inverse of lifetimes
        ThisLambda = sampleLambda(Data,EmpParam,S,...
            ThisLambda,Alpha,Beta,Alpha_Prop,N,T_IRF,T_intv,Sig_IRF);
         
    %Uniform profiles    
    elseif UnifRho == 1
        
        %Sampling profiles
        for nn = 1:3
            [ThisRho_IndP,AcceptGP]=sampleRho_Unif(Dtmp,EmpParam,S,ThisRho_IndP,AcceptGP);
        end
        for mm = 1:M
            ThisRho(mm).Rho = ThisRho_IndP(mm).Rho*ones(size(Cent_Confocal,1),1); 
        end         
        %Sampling assignment parameters
        S = sampleTag_Unif(Dtmp,EmpParam,ThisLambda,ThisRho_IndP,N);
        %Sampling inverse of lifetimes
        for nn = 1:3
            ThisLambda = sampleLambda_Unif(Dtmp,EmpParam,S,...
                ThisLambda,Alpha,Beta,Alpha_Prop,N,T_IRF,T_intv,Sig_IRF);
        end
         
    end
    
    if jj/PerSample == floor(jj/PerSample)
        Ind = Ind + 1;
        Chain(Ind).Rho_IndP = ThisRho_IndP;
        Chain(Ind).Rho = ThisRho;
        Chain(Ind).Lambda = ThisLambda;
        Chain(Ind).Xi = ThisXi;
    end
    
    %Saving every 50K chunk of the chain.
    if jj/50000 == floor(jj/50000)
        Num = (jj/50000-1);
        SInd = Num*floor(50000/PerSample)+1;
        tChain = Chain(SInd:end);
        save(sprintf('%s_%d',Str,Num),'tChain','S','-v7.3');
        sprintf('Part %d of the chain is saved\n',Num);
    end
      
end

end

