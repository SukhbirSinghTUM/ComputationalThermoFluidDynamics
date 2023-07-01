%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the parameters of the simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Geometry:
%
%    |---
%    |   ---
%    |      ----
%    |          |
% h1 |----------|  h2  <- symmetry axis
%    |          |
%    |      ----
%    |   ---
%    |---
%
%    |<--  l -->|


% Defin dimension of the trapezoidal domain
% h2 <= h1 !


shape = 'linear';  % 'linear' or 'quadratic'

h1 = 4;
hm = 4;            % only necessary for quatratic option 
h2 = 4;
l = 5;

% Number of degrees of freedom (number of nodes per length)
dimX = 50;
dimY = 50;


% Values for Dirichlet BC

TD.ambient = 298;
TD.htc = 50;

TD.flux = 3.5134e+07;


% Shape of the Cooling Fin

% h2 <= h1 !

switch shape
    
    case 'linear'

        formfunction = @(xnorm) (1-xnorm)*h1/2 + xnorm*h2/2;

    case 'quadratic'

        c1 = h2+2*h1/2-2*hm;
        c2 = 2*hm - 3*h1/2 - h2/2;
        c3 = h1/2;

        formfunction = @(xnorm) c1*xnorm.^2 +c2*xnorm + c3;
        
    case 'crazy'

        d1 = 3;
        d2 = 4;
        
        formfunction = @(xnorm) (1-xnorm)*h1/2 + xnorm*h2/2+ (sin(2*pi*d1*xnorm)).*(1-(1-1/d2)*xnorm);
        
    % Add other cases of form function
end



