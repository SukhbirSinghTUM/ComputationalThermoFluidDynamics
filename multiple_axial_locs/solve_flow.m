%% Solution for flow

%% Initialize variables
Init_flow

%% initialize spatial Matrix T
M_flow = zeros(dimY,dimX);

%% set up the mesh
[X_flow, Y_flow] = setUpMesh(M_flow, l, formfunction);

%% Solving for velocity
u = solveFVM_velocity(M_flow, X_flow, Y_flow);

% Temperature profile at the previous section in cooling channel
T_prev = 300*ones(dimY,dimX);

% Temperature profile at subsequent section 
T_new = solveFVM_temp(M_flow, X_flow, Y_flow, T_prev, u, Q);


%% Make some plots

figure(1)
surf(X_flow,Y_flow,u,'FaceColor','interp','EdgeColor','none')
ylim([Y_flow(end,1) Y_flow(1,1)])
view(0,90)
title('Velocity distribution')
colormap("turbo")
xlabel('x')
ylabel('y')
colorbar

figure(2)
surf(X_flow,Y_flow,T_new,'FaceColor','interp','EdgeColor','none')
ylim([Y_flow(end,1) Y_flow(1,1)])
view(0,90)
title('Temperature distribution')
colormap("turbo")
xlabel('x')
ylabel('y')
colorbar


figure(3)
contourf(X_flow,Y_flow,T_new,30)
ylim([Y_flow(end,1) Y_flow(1,1)])
colormap("turbo")
title('Temperature distribution')
xlabel('x')
ylabel('y')
colorbar
