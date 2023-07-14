%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the parameters of the flow part
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

h1_flow = 0.3; 
h2_flow = 0.3;
l_flow = 0.1;

% Number of degrees of freedom (number of nodes per length)
dimX = 50;
dimY = 50;

% Shape of the domain 
formfunction_flow = @(xnorm) (1-xnorm)*h1_flow/2 + xnorm*h2_flow/2;
