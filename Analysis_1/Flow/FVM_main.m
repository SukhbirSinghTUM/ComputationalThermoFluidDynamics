clc
clear 
close all

%% Initialize variables
InitFVM

%% initialize spatial Matrix T

M = zeros(dimY,dimX);

%% set up the mesh

[X, Y] = setUpMesh(M, l, formfunction);

%% Solving for velocity
u = solveFVM_velocity(M, X, Y);

% Temperature profile at the previous section in cooling channel
T_prev = 300*ones(dimY,dimX);

% Temperature profile at subsequent section 
T_new = solveFVM_temp(M, X, Y, T_prev, u);


%% Make some plots

figure(1)
surf(X,Y,u,'FaceColor','interp','EdgeColor','none')
ylim([Y(end,1) Y(1,1)])
view(0,90)
title('Velocity distribution')
colormap("turbo")
xlabel('x')
ylabel('y')
colorbar

figure(2)
surf(X,Y,T_new,'FaceColor','interp','EdgeColor','none')
ylim([Y(end,1) Y(1,1)])
view(0,90)
title('Temperature distribution')
colormap("turbo")
xlabel('x')
ylabel('y')
colorbar


figure(3)
contourf(X,Y,T_new,30)
ylim([Y(end,1) Y(1,1)])
colormap("turbo")
title('Temperature distribution')
xlabel('x')
ylabel('y')
colorbar
