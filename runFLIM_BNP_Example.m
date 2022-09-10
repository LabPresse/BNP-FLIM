%% Bayesian NonParametric (BNP) analysis of Fluorescence Lifetime Imaging (FLIM) Data 
%  This example generates data from a mixture of two lifetime components
%  with non-uniform profiles over a region of 10x10 pixels and 100 photons
%  per pixel. It then analyzes the generated data using the BNP-FLIM
%  technique and saves the results.
%
% Requirements and Setup:
%   1. MATLAB 2020 or higher versions
%   2. Statistics and Machine learning toolbox
%   3. BNP-FLIM software package
% 
% Results:
%   1. Chains of sampled parameters containing:
%       a) Chain of sampled inverse of lifetimes (Lambda)
%       b) Chain of sampled profiles at inducing points/pixel centers (Rho_IndP)
%       c) Chain of sampled profiles at test points (Rho)
%       c) Chain of sampled means of GP priors
%      Chain is splited into different chunks and saved independently.
%   2. Histogram of sampled lifetimes
%   3. Cross sections of the resulting profiles along the x-axis
%   

%% Populating EmpParams structure 
%EmpParam contains parameters of the experiment described in the following:
%These parameters are used in data generation and some of them are used in
%the algorithm.
%
% PixelSize: pixel size (mu m)
% Xi_X:      Confocal X center (mu m)
% Yi_X:      Confocal Y center (mu m)
% Zi_X:      Confocal Z center (mu m)
% Mu:        Array of excitation rates of differrent species (1/ns)
% Tau:       Array of lifetimes of the species (ns)
% Dp:        Average width of laser pulses (ns)
% Sig_IRF:   Sigma of IRF (ns)
% T_IRF:     Mean of IRF (ns)
% T:         Time interval between two consecurive laser pulses (ns)
% OmegaX:    HFWM of confocal PSF along X and Y axis (mu m)
% OmegaZ:    HFWM of confocal PSF along Z axis (mu m)   
% PSF:       The combination of detection and illumination PSFs (inline function)   
%

close all;clear;

M = 2; %Number of lifetimes components (species)
XPixNum = 10; %Number of pixels along the x-axis 
YPixNum = 10; %Number of pixels along the y-axis
EmpParam.PixelSize = 0.2; %Pixel size
EmpParam.Xi_X = [];
for pp = 1:XPixNum
    EmpParam.Xi_X = cat(1,EmpParam.Xi_X,EmpParam.PixelSize*(pp-0.5)*ones([YPixNum,1]));
end
EmpParam.Xi_Y = EmpParam.PixelSize*repmat((0.5:YPixNum-0.5)',[XPixNum,1]);
EmpParam.Xi_Z = zeros(XPixNum*YPixNum,1);
EmpParam.Mu = 0.001*ones(1,M);
%EmpParam.Tau = [0.8 3.6];
EmpParam.Tau = [1 15];
EmpParam.Dp = 0.1;
EmpParam.Sig_IRF = 0.66;
EmpParam.T_IRF = 10.4;
EmpParam.T = 15;
EmpParam.OmegaX = 1.4*EmpParam.PixelSize; 
EmpParam.OmegaZ = 3.5*EmpParam.PixelSize; 
EmpParam.PSF = @(X,Y,Xp,Yp,OmegaX,OmegaZ) OmegaZ*exp(-2*((X-Xp).^2/OmegaX^2 + ...
    (Y-Yp).^2/OmegaX^2));

%% Generate Data
% The Data structure contains the following parameters.
% Dt:         Photon arrival times (microtimes) (ns)
% DtTot:      Time intervals from excitation to detection of photons
% W:          Binary parameter determining empty and non-empty pulses
% S:          Assignment parameter associating photons to species (This  
%             parameter is not used in the analysis. The algorithm learns it)
% X_Confocal: X-locations of the confocal/pixel centers
% Y_Confocal: Y-locations of the confocal/pixel centers
% Z_Confocal: Z-locations of the confocal/pixel centers

addpath('BNP_Functions')

% lifetime profiles 
RhoIn(1).Rho = @(X,Y,Z) 1000+1000*erfc((X-XPixNum*EmpParam.PixelSize/2)*2);
RhoIn(2).Rho = @(X,Y,Z) 1000+1000*erfc(-(X-XPixNum*EmpParam.PixelSize/2)*2);

D = EmpParam.PixelSize; %Pixel size
Np = 200; %Number of photons per pixel

%Generating data
Data = genFLIM(EmpParam,RhoIn,D,Np);

%% BNP structures contains parameters used in the algorithm
% M:            Number of lifetime components (species)
% Alpha:        Shape parameter of gamma prior on inverse of lifetimes
% Beta:         Scale parameter of gamma prior on inverse of lifetimes
% Alpha_Prop:   Parameter of proposal dist. of inverse of lifetimes
% NJump:        Number of samples from the posterior 
% D:            Grid size of the test points (resolution of profiles) (mu m)
% T:            Parameter of GP prior
% L:            Parameter of GP prior (range of correlation) (mu m)
% N:            Number of preceeding laser pulses to be considered
% Sig_GP:       Average jump size in proposing new profiles
% Sig_Xi:       Average jump size in proposing new Xi
% Sig_Prior_Xi: Sigma of Gaussian prior on Xi 
% DEBUG:        Show the current profiles (1) (Default=0)
% PerSample:    Adding one element in the chain per PerSample samples (Default=50) 
% Exten:        Extending profiles at the edges during inference (Default=2 pixels)
%               

BNP.M = M;
BNP.Alpha = 1; 
BNP.Beta = 25; 
BNP.Alpha_Prop = 5000; 
BNP.NJump = 350000; %More samples may need to be taken for different problems
BNP.NJump = 500000; 
BNP.D = EmpParam.PixelSize/2; 
BNP.T = 1; 
BNP.L = 0.8;
BNP.N = 4; 
BNP.Sig_GP = 0.007;
BNP.Sig_Xi = 0.5;
BNP.Sig_Prior_Xi = 50;
BNP.DEBUG = 0;

%% Making Inference
%This part is time consuming and will probably take ~3 hours 

Data = Data(:);
EmpParam.Mu = ones(1,BNP.M); 
EmpParam.Dp = 1;

tic;
Chain=runBNPs_FLIM(Data,EmpParam,BNP);
T = toc;
fprintf('It took %d s to analyze this data set.\n',T)

%% Generating Results using the last 1/8th of the output chain
%This section generates a plot containing two subplots. The first subplot
%depicts a cross-section of the simulated profiles along the x-axis. The
%second subplot shows a histogram of the lifetime samples.
%

Tau = [];
Rho1 = [];
Rho2 = [];
NPt = sqrt(length(Chain(end).Rho(1).Rho));
for ii = 5001:6000
    Tau = cat(1,Tau,1./Chain(ii).Lambda);
    tRho1 = reshape(Chain(ii).Rho(1).Rho,[NPt NPt]);
    tRho2 = reshape(Chain(ii).Rho(2).Rho,[NPt NPt]);
    Rho1 = cat(1,Rho1,tRho1(5,5:24));
    Rho2 = cat(1,Rho2,tRho2(5,5:24));
end

figure('unit','inches','Position',[5 7 7.5 3.5]);subplot(1,2,2)
histogram(Tau,'BinWidth',0.1,'normalization','pdf','orientation','horizontal')
hold;plot(0:8,EmpParam.Tau(1)*ones(1,9),'--k','linewidth',1.6);
X = 0:0.1:8;Y = gampdf(1./X,BNP.Alpha,BNP.Beta);
plot(Y,X,'linewidth',1.6)
plot(0:8,EmpParam.Tau(2)*ones(1,9),'--k','linewidth',1.6)
ylim([0 max(EmpParam.Tau)+2]);xlim([0 6]);set(gca,'Ytick',EmpParam.Tau,'FontSize',12);set(gca,'Xtick',[0 6])
ylabel(sprintf('\\tau (ns)'),'FontSize',12)
xlabel('pdf','FontSize',12);box off
legend({'lifetime samples','Ground truth','Prior'})
legend boxoff;

Xg = BNP.D/2:BNP.D:EmpParam.PixelSize*XPixNum;
Rho1True = RhoIn(1).Rho(Xg,Xg,0)*0.1;
Rho2True = RhoIn(2).Rho(Xg,Xg,0)*0.1;
Xfill = [Xg,fliplr(Xg)];
Between1 = [quantile(Rho1,0.0025),fliplr(quantile(Rho1,0.9975))];
subplot(1,2,1);fill(Xfill,Between1(1,:,:)*EmpParam.Dp,[0 0.4470 0.9410],'FaceAlpha',0.4,'LineStyle','none')
Between2 = [quantile(Rho2,0.0025),fliplr(quantile(Rho2,0.9975))];
hold;plot(Xg,EmpParam.Dp*0.001*Rho1True,'--k','linewidth',1.8);
plot(Xg,EmpParam.Dp*median(Rho1),'m','linewidth',1.5);
plot(Xg,EmpParam.Dp*0.001*Rho2True,'k--','linewidth',1.8)
fill(Xfill,Between2(1,:,:)*EmpParam.Dp,[0 0.4470 0.9410],'FaceAlpha',0.4,'LineStyle','none')
plot(Xg,EmpParam.Dp*median(Rho2),'m','linewidth',1.5);plot(Xg,EmpParam.Dp*median(Rho1),'m','linewidth',1.5);
ylim([0 0.53]);xlim([0 2]);set(gca,'Ytick',[0 0.1 0.3],'FontSize',12);set(gca,'Xtick',[0 1 2])
Ay=ylabel(sprintf('\\mu\\rho (%sm^{-3})',char(956)),'FontSize',12);box off;
Ay.Position(1) = -0.42;
xlabel(sprintf('x (%sm)',char(956)),'FontSize',12)
box off;
legend({'95% confidence interval','Ground truth','Median'})
legend boxoff;
