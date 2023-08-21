% This function calculates the heat flux going to the solid side from hot
% gas side at a given cross section of nozzle. This particular function is
% for LO2 and RP1 rockets. 
% The inputs are:
%           1. A:   Area of cross section being studied
%           2. A_t: Throat area of the nozzle being considered
%           3. Tw:  The vector of wall temperatures on solid side


function Q = heat_flux(A, A_t, Tw)

% Stagnation Conditions For LO2 and RP-1
Pr = 1;                     % Prandtl number         
Cp = 1.881;                 % J/kgK 
mu = 1e-3;                % m2/sec

D_t = sqrt(4*A_t/pi);
P0 = 1e7;                  % Stagnation Pressure in chamber: 10Mpa
g = 9.8;                    % Acc due to gravity
k = 1.2;                    % Ratio of specific heat capacities of exhaust gases

T0 = 3273;                                                   % Chamber temperature: 3000 degC
R = 0.3e03;                                                  % Gas constant: 0.3kJ/kgK
c_star = sqrt((R*T0/k)*((0.5*(k+1))^((k+1)/(k-1))));         % Characteristic Velocity

r = 5e-03;                  % Radius of curvature at throat: R = 5mm
 
% Mach number from given area of cross section using mach-area relation
fcn = @(Ma)((1/Ma^2)*((((2/(k+1))*(1 + 0.5*(k-1)*Ma^2))^((k+1)/(k-1))) - (A/A_t)^2));
Ma = fzero(fcn,1);

% Correction factor: sigma
K1 = ((0.5.*(Tw./T0).*(1 + 0.5*(k-1)*(Ma^2))) + 0.5).^(-0.68);
K2 = (1+0.5*(k-1)*(Ma^2))^(-0.12);
sigma = K1.*K2;

% Heat Transfer Coefficient 
C1 = (mu^0.2)*Cp/(Pr^0.6);
C2 = (P0*g)/c_star;
htc = (0.026/(D_t^0.2))*C1*(C2^0.8)*((D_t/r)^0.1)*((A_t/A)^0.9).*sigma;

% Calculating Gas Temperature from adiabatic relations
Tg = T0*(1+0.5*(k-1)*(Ma^2))^-1;

% Heat flux
Q = 5.*htc.*(Tg - Tw);

end 