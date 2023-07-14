clc
clear 
close all


%% Initialize variables

InitFVM

% Physical properties of nozzle
D_t = 0.15;                 % Throat Diameter: 15cm 
A_t = pi*(D_t^2)/4;         % Throat Area

%% initialize spatial Matrix T

M = zeros(dimY,dimX);

%% set up the mesh

[X, Y] = setUpMesh(M, l, formfunction);


%% Iterating 

% Initial guess for Tw wall temperature (ignoring corners)
Tw = 1000.*ones(1,size(M,2));


eps = 1;                             % For stopping criteria
idx = 1;                             % index number for the coming loop

max_Tw_h = zeros(1,100);             % Storing the max value of temperature at wall on hot side       
max_Tw_c = zeros(1,100);             % Storing the max value of temperature at wall on solid side 
fprintf('Iteration      |       Residual\n');

for i = 1:100

    Tw_h = Tw;                                   % For the purpose of comparing the hotter and colder wall sides
    Q = 1e-4.*heat_flux(10*A_t, A_t, Tw_h);      % Heat flux at every southern node except corners
    
    T = solveFVM(M, X, Y, TD, Q);
    T = reshape(T,dimY,dimX);
    
    Tw_c = T(size(M,1),:);                      % Extracting the temp at southern nodes in solid domain
    Tw = real(Tw_c);                            % Keeping only the real parts of solution

    max_Tw_h(idx) = max(Tw_h);                  % Storing the max value of temp at hot side in max_Tw_h vector
    max_Tw_c(idx) = max(Tw_c);                  % Storing the max value of temp at solid side in max_Tw_c vector

    fprintf('%d              |       %f \n',idx,norm(Tw_h - Tw_c));

    % Checking if the solution has converged or not
    if norm(Tw_h - Tw_c)<eps                
        break
    elseif i == 100
        fprintf("Max iterations reached");
    end 

    idx = idx + 1;

end 


%% Make some plots

% Plotting the residuals vs no. of iterations
figure(1)
plot(1:idx,max_Tw_h(1:idx),'r-o',1:idx,max_Tw_c(1:idx),'b-o')
xlim([1 idx])
xlabel('Iteration','FontSize',14)
ylabel('Wall Temperature','FontSize',14)
legend('Hot gas side','Solid side','FontSize',14)

figure(2)
surf(X,Y,T,'FaceColor','interp','EdgeColor','none')
view(0,90)
title('Temperature distribution')
colormap("turbo")
xlabel('x')
ylabel('y')
colorbar

