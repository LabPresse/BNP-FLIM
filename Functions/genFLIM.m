function Data = genFLIM(EmpParam,RhoIn,D,Np,Flag)
%genFLIM generates FLIM data given species concentrations.
%
%INPUT:
%   EmpParam: Structure containing parameters required to simulate the data
%       Xi_X:   Confocal X center (pixel)
%       Yi_X:   Confocal Y center (pixel)
%       Zi_X:   Confocal Z center (pixel)
%       Mu:     Array of excitation rate of the species (1/ns)
%       Tau:    Array of lifetimes of the species (ns
%       DeltaP: Pulse duration (ns)
%       Sig_IRF:Sigma of arrival times error (ns)
%       T_IRF:  Mean of arrival times errors (ns)
%       OmegaX: HFWM of confocal PSF along X and Y axis (pixel)
%       OmegaZ: HFWM of confocal PSF along Z axis (pixel)   
%       PSF:    The combination of detection and illumination PSFs (inline function)       
%   RhoIn:  Structure containing inline functions of concentrations of 
%           molecular species.
%   Np:     Number of photons per confocal spot
%
%OUTPUT:
%   Data:   Structure array where an element contains photon arrival times 
%           for corresponding confocal spot. (ns)
%
%Created by:
%   Mohamadreza Fazel (Presse Lab, 2020)
%

if nargin < 5
   Flag = 'Photon'; 
end

if ~strcmp(Flag,'Photon') && ~strcmp(Flag,'Pulse')
    error('The Flag must be either Photon or Pulse.');
end

Exten = 2;
XPixNum = max(EmpParam.Xi_X)/EmpParam.PixelSize + 0.5;
YPixNum = max(EmpParam.Xi_Y)/EmpParam.PixelSize + 0.5;
MinRangeX = -(Exten*EmpParam.PixelSize-D/2)*[1, 1]; 
MaxRangeX(1) = (XPixNum*EmpParam.PixelSize+Exten*EmpParam.PixelSize);
MaxRangeX(2) = (YPixNum*EmpParam.PixelSize+Exten*EmpParam.PixelSize);
[Xg,Yg,Zg] = meshgrid(MinRangeX(1):D:MaxRangeX(1),MinRangeX(2):D:MaxRangeX(2),0);
TestP = [Xg(:),Yg(:),Zg(:)];
M = length(RhoIn);
for mm = 1:M
    RhoInT(mm).Rho = zeros(length(Xg(:)),1);
end

for mm = 1:M
    Rhotmp = RhoIn(mm).Rho;
    if isscalar(Rhotmp)
        Rho_init = Rhotmp(Xg(:),Yg(:),Zg(:));
        RhoInT(mm).Rho = Rho_init;
    else
        RhoInT(mm).Rho = RhoIn(mm).Rho(:);
    end
end

if M ~= length(EmpParam.Tau) || M ~= length(EmpParam.Mu)
   error('Numbers of input concentration profiles, lifetimes and mu must be the same.'); 
end

% Rhotmp = RhoIn(3).Rho;
% Rho3_init = Rhotmp(Xg,Yg,Zg);
% RhoInT(3).Rho = Rho3_init(:);

LXi = length(EmpParam.Xi_X);
Data(LXi).Dt = [];
Data(LXi).DtTot = [];
Data(LXi).W = [];
Data(LXi).S = [];
Data(LXi).X_Confocal = [];
Data(LXi).Y_Confocal = [];
Data(LXi).Z_Confocal = [];
for ii = 1:LXi
    Data(ii).X_Confocal = EmpParam.Xi_X(ii);
    Data(ii).Y_Confocal = EmpParam.Xi_Y(ii);
    Data(ii).Z_Confocal = EmpParam.Xi_Z(ii);
end

Pi_A = zeros(LXi,M);
P0A = zeros(1,M);
%Pi_B = zeros(LXi,1); 
%Pi_AB = zeros(LXi,1);

for ii = 1:LXi
    for mm = 1:M
        IntA = discreteInt(TestP,RhoInT(mm).Rho,EmpParam,D,ii);
        P0A(mm) = exp(-EmpParam.Mu(mm)*EmpParam.Dp*IntA);
    end
    
    for mm = 1:M
        Pi_A(ii,mm) = (1-P0A(mm))*prod(P0A)/P0A(mm);
    end
end
T = EmpParam.T;
Pi_Sum = cumsum(Pi_A,2);
for jj = 1:LXi
    ii = 0;
    NPulse = 0;
    This_P = Pi_Sum(jj,:);
    while ii < Np
        if strcmp(Flag,'Pulse')
            ii = ii + 1;
        end
        NPulse = NPulse + 1;
        Rnd = rand();
        Ind=find(This_P-Rnd>0,1);
        if ~isempty(Ind)
            if strcmp(Flag,'Photon')
                ii = ii + 1;
            end
            Data(jj).W = cat(1,Data(jj).W,1);
            Ttmp = exprnd(EmpParam.Tau(Ind))+normrnd(EmpParam.T_IRF,EmpParam.Sig_IRF);
            Data(jj).Dt = cat(1,Data(jj).Dt,Ttmp-T*floor(Ttmp/T));
            Data(jj).DtTot = cat(1,Data(jj).DtTot,Ttmp);
            Data(jj).S = cat(1,Data(jj).S,Ind);
        else
            Data(jj).W = cat(1,Data(jj).W,0);
        end
    end
    
    if sum(Data(jj).W) == 0
       error('There is a problem in generating the data');
    end
    
end
end