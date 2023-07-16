function [T, Q_flow_solid2, Q_flow_SOLID] = solveFVM_temp(M, X, Y, T_prev, u, Q, Tw_flow_solid2, Tw_SOLID_flow)


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
        [A(index(i,j), :), B] =  stamp_temp(i, j, X, Y, B, T_prev, u, Q, Tw_flow_solid2, Tw_SOLID_flow);     
    end
end


% solve the linear system
A = sparse(A);

T = A\B;
T = reshape(T,size(M,1),size(M,2));

k_fluid = 0.12;              % Thermal conductivity of nickel alloys = 300W/mK
Q_flow_solid2 = abs(k_fluid.*(T(1,:)- T(2,:))./(Y(1,:)-Y(2,:)));
Q_flow_SOLID = abs(k_fluid.*(T(:,1)- T(:,2))./(X(:,1)-X(:,2)));







