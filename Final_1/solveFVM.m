function [T, Q_out1, Q_out2] = solveFVM(solid_type, M, X, Y, Tw_solid_flow_s, Tw_west ,Q_south, Q_east )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File solveFVM.m
%
% This routine set up the linear system and solve it
%
% input
% M         Spatial Matrix M
% X         Matrix x coordinates 
% Y         Matrix y coordinates
% boundary  String vector. Boundary types.
% TD        Temperature for each boundary (if Dirichlet)
% alpha     convective heat transfer coefficient
% Tinf      Temperature of the surrouding fluid 
%
% output
% T         Temperature field

% Index maps the node position to the correct linear equation

index = @(ii, jj) ii + (jj-1) * size(M, 1);


% B is the right-hand side of the linear system
B = zeros(size(M,1)*size(M,2),1);

% set up the system matrix A
A = zeros(size(M, 1)*size(M, 2),size(M, 1)*size(M, 2));

for i = 1:size(M, 1)
    for j = 1:size(M, 2)
        % Fill the system matrix and the right-hand side for node (i,j)
        [A(index(i,j), :), B] =  stamp_mex(i, j, X, Y, B, solid_type, Tw_solid_flow_s, Tw_west , Q_south, Q_east);     
    end
end


% solve the linear system
A = sparse(A);

T = A\B;
T= reshape(T,size(M,1),size(M,2));

if strcmp(solid_type,"solid1")
    k_solid = 0.1;              % Thermal conductivity of nickel alloys = 300W/mK
    Q_out1 = -(k_solid.*(T(1,:)- T(2,:))./(Y(1,:)-Y(2,:)));
    Q_out2 = -(k_solid.*(T(:,1)- T(:,2))./(X(:,1)-X(:,2)));
end 

if strcmp(solid_type,"solid2")
    k_solid = 0.1;              % Thermal conductivity of nickel alloys = 300W/mK
    Q_out1 = [];
    Q_out2 = -(k_solid.*(T(1,1)- T(1,2))./(X(:,1)-X(:,2)));
end 

if strcmp(solid_type,"SOLID")
    Q_out1 = [];
    Q_out2 = [];
end 






