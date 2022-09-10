function Integ = discreteInt(TestP,Rho,EmpParam,D,ii,IDTestP)
%discreteInt finds integral of Rho using the given discrete values.
%
%INPUT: 
%   TestP:    Set of discrete test points
%   Rho:      Concentrations
%   EmpParam: Structure containing parameters of the experiment
%   D:        Separation between the test points
%   ii:       Pixel number
%   IDTestP:  Set of flags determining the test points no further than
%             3sigma from the pixel center.
%
%OUTPUT:
%   Integ:    The calculated sum (Integral)
%
%Created by:
%   Mohamadreza Fazel (Presse lab 2020)
%

Xp = EmpParam.Xi_X(ii);
Yp = EmpParam.Xi_Y(ii);
Zp = EmpParam.Xi_Z(ii);
if nargin < 6 || isempty(IDTestP) 
    Dis = pdist2(TestP,[Xp,Yp,Zp]);
    ID = Dis<3*EmpParam.OmegaX;
else
   ID = IDTestP(:,ii)~=0;  
end
TestP_tmp = TestP(ID,:);

Integ=sum(Rho(ID).*EmpParam.PSF(TestP_tmp(:,1),TestP_tmp(:,2),Xp,Yp,...
    EmpParam.OmegaX,EmpParam.OmegaZ))*D^2;

end

