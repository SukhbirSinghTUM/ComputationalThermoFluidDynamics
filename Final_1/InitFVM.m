%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the parameters of the Solid-1 part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Geometry:
%
%    |----------|
%    |          |   
% h1 |          | h2   
%    |          |
%    |----------|  
%    |<--  l -->|

shape = 'linear';  % 'linear' or 'quadratic'

h1_solid1 = 0.15;
h2_solid1 = 0.15;
h1_SOLID = 2*h1_solid1+h1_flow;
h2_SOLID = 2*h2_solid1+h2_flow;
l_solid1 = 0.1;


% Number of degrees of freedom (number of nodes per length)
dimX = 50;
dimY = 50;

% Shape of the domain
formfunction_solid1 = @(xnorm) (1-xnorm)*h1_solid1/2 + xnorm*h2_solid1/2;
formfunction_SOLID = @(xnorm) (1-xnorm)*h1_SOLID/2 + xnorm*h2_SOLID/2;

