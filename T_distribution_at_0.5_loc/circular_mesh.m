% This function gives out the mesh by taking in the following inputs:
%   1. M        : spatial matrix containing the dimX and dimY for discretization
%   2. r        : inner radius of the nozzle
%   3. h        : height upto which the mesh should go (radial height)
%   4. theta    : angle defining the thickness of mesh

%       \                         /
%        \                       /
%         \____________________ /
%          \                   /   ^  
%           \                 /    |  
%            \               /     | h    
%             \ ___________ /      |
%              \           /    ^
%               \         /     |
%                \       /      |
%                 \theta/       | r
%                  \   /        |
%                   \ /         |    



function [X , Y] = circular_mesh(M, r, h, theta)

dimR = size(M,1);
dim_theta = size(M,2);

R1 = r ;         % inner radius 
R2 = r + h ;     % outer radius
nR = linspace(R2, R1, dimR) ;
nT = linspace(0,theta,dim_theta);
[R, T] = meshgrid(nR,nT) ;

X = R.*sind(T); 
Y = R.*cosd(T);


X = X';
Y= Y';

end