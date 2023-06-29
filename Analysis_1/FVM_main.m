clc
clear 
close all


%% Initialize variables

InitFVM

%% initialize spatial Matrix T

M = zeros(dimY,dimX);

%% set up the mesh

[X, Y] = setUpMesh(M, l, formfunction);

%% Iterating 

T = solveFVM(M, X, Y, TD);
T = reshape(T,dimY,dimX);


%% Make some plots
figure(1)
surf(X,Y,T,'FaceColor','interp')
title('Temperature distribution')
colormap("turbo")
xlabel('x')
ylabel('y')
colorbar

figure(2)
contourf(X,Y,T,30)
colormap("turbo")
title('Temperature distribution')
xlabel('x')
ylabel('y')
colorbar
