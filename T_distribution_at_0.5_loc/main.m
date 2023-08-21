clc
clear 
close all

%% Initialize variables

% discretization of solid1, flow and solid2 have 20X20 nodes each
dimX = 20;
dimY = 20;

% Physical properties of nozzle
D_t = 0.05;                 % Throat Diameter: 5cm 
A_t = pi*(D_t^2)/4;         % Throat Area
D_e = 0.8;                  % Exhaust Diameter: 80cm 
A_e = pi*(D_e^2)/4;         % Exhaust Area

% Thickness of sub-domains
h1 = 0.01;                  % thickness of solid1 subdomain
h2 = 0.01;                  % thickness of flow subdomain
h3 = 0.01;                  % thickness of solid2 subdomain

%% initialize spatial Matrix T
M_solid1 = zeros(dimY,dimX);       
M_flow = zeros(dimY,dimX);
M_solid2 = zeros(dimY,dimX);

% spatial matrix of SOLID (right half) will have two lesser rows because of
% the interfaces with right half. 
M_SOLID = zeros(size(M_solid1,1)+size(M_flow,1)+size(M_solid2,1)-2,dimX);


%% Temperature Inlet for flow
% The coolant enters from the nozzle exit at 300K or room temperature
T_prev = 300.*ones(dimY,dimX);

%% Initial guesses starting from axial location = 1

% Initial guess for hot-gas wall temperature 
Tw_gas = 0.*ones(1,size(M_solid1,2)+size(M_SOLID,2));           % Temp of interface joining gas and solid1+SOLID
Tw_gas_solid1 = Tw_gas(1,end - size(M_solid1,2) + 1:end);       % Temp of interface joining gas and solid1
Tw_gas_SOLID = Tw_gas(1, 1:size(M_SOLID,2));                    % Temp of interface joining gas and SOLID

% Initial guess for solid1+flow wall temperature
Tw_solid1_flow = 0.*ones(1,size(M_flow,2));

% Initial guess for flow+solid2 wall temperature
Tw_flow_solid2 = 0.*ones(1,size(M_flow,2));

% Initial guess for SOLID+right wall temperature
Tw_SOLID = 0.*ones(1,size(M_SOLID,1))';

%% Loop for axial locations
axial_loc = 1:-0.1:0.1;                     % location of axial points being studied       
max_iter = 400;                             % maximum number of iterations allowed 

% looping over all the axial locations (necessary to reach at the middle
% axial_loc = 0.5), because the temp in cooling channel develops from
% axial_loc = 1 to axial_loc =0

for loc = axial_loc
    
    % Here we find out the area of nozzle at diff axial locations
    % This also gives the pressure gradient to be considered at the location
    [Area_diff, P_dif] = complete_nozzle(D_t, loc, D_e);
    
    %% set up the mesh
    R = 0.5*sqrt((4/pi)*Area_diff);             % Radius of nozzle at axial loc
    
    % Setting up the circular mesh 
    [X_solid1, Y_solid1] = circular_mesh(M_solid1, R ,h1, 3);
    [X_flow, Y_flow] = circular_mesh(M_flow, R+h1 ,h2, 3);
    [X_solid2, Y_solid2] = circular_mesh(M_solid2, R+h1+h2 ,h3, 3);
    
    % Here we combine the mesh from right half of domain and combine them
    % to form the "SOLID" mesh
    [X_SOLID1, Y_SOLID1] = circular_mesh(M_solid1, R ,h1, -3);
    [X_SOLID2, Y_SOLID2] = circular_mesh(M_flow, R+h1 ,h2, -3);
    [X_SOLID3, Y_SOLID3] = circular_mesh(M_solid2, R+h1+h2 ,h3, -3);
    
    X_SOLID = [X_SOLID3; X_SOLID2; X_SOLID1];
    Y_SOLID = [Y_SOLID3; Y_SOLID2; Y_SOLID1];
    X_SOLID([size(M_solid2,1), size(M_solid2,1)+size(M_flow,1)],:) = [];
    Y_SOLID([size(M_solid2,1), size(M_solid2,1)+size(M_flow,1)],:) = [];
    X_SOLID = flip(X_SOLID,2);
    Y_SOLID = flip(Y_SOLID,2);
    
    %% Solving velocity part of flow
    % Here the axial velocity at every node of fluid mesh is found out
    u = solveFVM_velocity(M_flow, X_flow, Y_flow, P_dif);
    
    %% Iterating at current section

    eps = 0.01;                           % For stopping criteria
    
    % These error vectors store the difference in two temperatures on
    % either side of the interface. example: error_gas_solid1 stores the
    % difference between temp. of wall on gas side and solid1 side 
    error_gas_solid1 = zeros(1,100);
    error_solid1_flow = zeros(1,100);
    error_flow_solid2 =   zeros(1,100);
    error_left_right =  zeros(1,100);
    
    fprintf('__________________________________________________________________________\n');
    fprintf('                                Residuals                                 \n');
    fprintf('__________________________________________________________________________\n\n');
    fprintf('Iteration\t|\tAxial Location\t|\tGas-Solid1\t|\tSolid1_Flow\t|\tFlow-Solid2\t|\tSOLID-right\n');
    
    
    residuals = [];
    for iter = 1:max_iter
    
        % Solving the Solid1 first, assuming some wall temperature on gas side and SOLID side
        solid_type = "solid1";
        Tw_gas_solid1_h = Tw_gas_solid1;                                   % Temperature on hot gas side equated to initial guess of wall temperature
        Tw_solid1_flow_s = Tw_solid1_flow;                                 % Temperature on flow side equated to initial guess of wall temperature
        Tw_solid1_SOLID = Tw_SOLID(end-dimY+1:end);                        % Temperature on SOLID side equated to initial guess of wall temperature
        Q_gas_solid1 = heat_flux(Area_diff, A_t, Tw_gas_solid1_h);         % Heat flux coming from hot-gas at every southern node of solid1 
        
        %Solving for temp profile in solid1
        [T_solid1, Q_solid1_flow, Q_solid1_SOLID] = solveFVM(solid_type, M_solid1, X_solid1, Y_solid1, Tw_solid1_flow_s, Tw_solid1_SOLID, Q_gas_solid1, [] );
        T_solid1 = real(T_solid1);                                         % ignoring imaginary components 
        Tw_gas_solid1_s = T_solid1(end,:);                                 % Extracting the temp at southern nodes in solid1 domain

        % updating the interface temperature of gas+solid1 
        Tw_gas_solid1 = Tw_gas_solid1_s;       
    
        %storing the error
        error_gas_solid1(1,iter) = max(abs(Tw_gas_solid1_h - Tw_gas_solid1_s));
        
        % Solving the flow part assuming some wall temperature of solid1+flow
        Tw_flow_solid2_f = Tw_flow_solid2;
        [T_flow, Q_flow_solid2, Q_flow_SOLID] = solveFVM_temp(M_flow, X_flow, Y_flow, T_prev, u, Q_solid1_flow, Tw_flow_solid2_f, Tw_SOLID(dimY + 1:2*dimY));
        T_flow = real(T_flow);
        T_prev = T_flow;
        Tw_solid1_flow_f = T_flow(end,:);

        % updating the interface temperture of solid1+flow
        Tw_solid1_flow = Tw_solid1_flow_f;
    
        error_solid1_flow(1,iter) = max(abs(Tw_solid1_flow_s - Tw_solid1_flow_f));
    
        % Solid type Solid2
        solid_type = "solid2";
        Tw_solid2_SOLID = Tw_SOLID(1:dimY);
        [T_solid2, ~, Q_solid2_SOLID] = solveFVM(solid_type, M_solid2, X_solid2, Y_solid2, [] ,Tw_solid2_SOLID , Q_flow_solid2,[]);
        T_solid2 = real(T_solid2);
    
        Tw_flow_solid2_s = T_solid2(end,:);

        % updating the interface temperture of flow+solid2
        Tw_flow_solid2 = Tw_flow_solid2_s; 
    
        error_flow_solid2(1,iter) = max(abs(Tw_flow_solid2_s - Tw_flow_solid2_f));
    
        % Solving the SOLID (left half) part
        Tw_gas_SOLID_h = Tw_gas_SOLID;                                      
        Q_gas_SOLID = heat_flux(Area_diff, A_t, Tw_gas_SOLID_h); 
        Q_east_SOLID = [Q_solid2_SOLID; Q_flow_SOLID; Q_solid1_SOLID];
        Q_east_SOLID([size(M_solid2,1), size(M_solid2,1)+size(M_flow,1)],:) = [];
    
        solid_type = "SOLID";
        Tw_SOLID_left = Tw_SOLID;
        [T_SOLID, ~, ~] = solveFVM(solid_type, M_SOLID, X_SOLID, Y_SOLID, [] ,[] , Q_gas_SOLID, Q_east_SOLID);
        T_SOLID = real(T_SOLID);
    
        Tw_SOLID = T_SOLID(:,end);

        % updating the interface temperture of gas+SOLID
        Tw_gas_SOLID = T_SOLID(end,:);
    
        error_left_right(1,iter) = max(abs(Tw_SOLID_left - T_SOLID(:,end)));
    
    
        %printing the residuals part
        iterationStr = sprintf('%-9d', iter);
        residual1 = sprintf('%.6f', max(abs(Tw_gas_solid1_h - Tw_gas_solid1_s)));
        residual2 = sprintf('%.6f', max(abs(Tw_solid1_flow_s - Tw_solid1_flow_f)));
        residual3 = sprintf('%.6f', max(abs(Tw_flow_solid2_s - Tw_flow_solid2_f)));
        residual4 = sprintf('%.6f', max(abs(Tw_SOLID - Tw_SOLID_left)));
        axial_location = sprintf('%.2f', loc);
    
        fprintf('%s\t|\t%s\t\t\t|\t%s\t|\t%s\t|\t%s\t|\t%s\n', iterationStr, axial_location, residual1, residual2, residual3,residual4);
    
    
        % Checking if the solution has converged or not
        if max(abs(Tw_gas_solid1_h - Tw_gas_solid1_s))<eps &&...
                max(abs(Tw_solid1_flow_s - Tw_solid1_flow_f)) <eps &&...
                max(abs(Tw_flow_solid2_s - Tw_flow_solid2_f)) <eps &&...
                max(abs(Tw_SOLID - Tw_SOLID_left)) <eps
            fprintf("Solution has converged.  \n")
            break
        elseif iter == max_iter
            fprintf("Max iterations reached.  \n");
        end 
    
    end 
    
    if loc == 0.5
        f1 = figure(1);
        f1.WindowState = 'maximized';
        hold on;
        plot(error_gas_solid1(1:iter),'LineWidth',2,'DisplayName','Gas+Solid1');
        plot(error_solid1_flow(1:iter),'LineWidth',2,'DisplayName','Solid1+Flow');
        plot(error_flow_solid2(1:iter),'LineWidth',2,'DisplayName','Flow+Solid2');
        plot(error_left_right(1:iter),'LineWidth',2,'DisplayName','SOLID+Right');
        xlim([1 iter+1])
        title('Residuals')
        xlabel('Iteration','FontSize',14)
        ylabel('Residuals','FontSize',14)
        legend('FontSize',14,'Location','northeast')
        drawnow
        
              
        %% Plot for the temperature distribution
        
        X_solid2(end,:)=[];
        Y_solid2(end,:)=[];
        T_solid2(end,:)=[];
        Y_flow(end,:)=[];
        X_flow(end,:)=[];
        T_flow(end,:)=[];
        
        Y_right = [Y_solid2; Y_flow; Y_solid1];
        X_right = [X_solid2; X_flow; X_solid1];
        T_right = [T_solid2; T_flow; T_solid1];
        Y_all = [Y_SOLID Y_right];
        X_all = [X_SOLID X_right];
        T_all = [T_SOLID T_right];
        drawnow
        
       
        f2 = figure(2);
        f2.WindowState = 'maximized';
        surf(X_all,Y_all,T_all,'FaceColor','interp','EdgeColor','none');
        view(0,90);
        ylim([min(Y_all,[],'all') max(Y_all,[],'all')]);
        title("Temperature distribution at axial location = "+ axial_location,'FontSize',15);
        colormap("turbo");
        xlabel('x','FontSize',15);
        ylabel('y','FontSize',15);
        colorbar;
        hold on
        plot3( [X_SOLID(1,1:end) (X_SOLID(1:end,end))' flip(X_SOLID(end,1:end)) flip((X_SOLID(1:end,1))')],...
            [Y_SOLID(1,1:end) (Y_SOLID(1:end,end))' flip(Y_SOLID(end,1:end)) flip((Y_SOLID(1:end,1))')],...
            5000.*ones(1,2*(size(X_SOLID,1)+size(X_SOLID,2))),'color','r','linewidth',1.5)
        plot3( [X_solid1(1,1:end) (X_solid1(1:end,end))' flip(X_solid1(end,1:end)) flip((X_solid1(1:end,1))')],...
            [Y_solid1(1,1:end) (Y_solid1(1:end,end))' flip(Y_solid1(end,1:end)) flip((Y_solid1(1:end,1))')],...
            5000.*ones(1,2*(size(X_solid1,1)+size(X_solid1,2))),'color','r','linewidth',1.5)
        plot3( [X_flow(1,1:end) (X_flow(1:end,end))' X_solid1(1,end)],...
            [Y_flow(1,1:end) (Y_flow(1:end,end))' Y_solid1(1,end)],...
            5000.*ones(1,(size(X_flow,1)+size(X_flow,2))+1),'color','r','linewidth',1.5)
        plot3( [X_solid2(1,1:end) (X_solid2(1:end,end))' X_flow(1,end)],...
            [Y_solid2(1,1:end) (Y_solid2(1:end,end))' Y_flow(1,end)],...
            5000.*ones(1,(size(X_solid2,1)+size(X_solid2,2))+1),'color','r','linewidth',1.5)
        drawnow

        u(end,:)=[];
        f3 = figure(3);
        f3.WindowState = 'maximized';
        surf(X_flow,Y_flow,u,'FaceColor','interp','EdgeColor','none');
        view(0,90);
        ylim([min(Y_flow,[],'all') max(Y_flow,[],'all')]);
        title("Flow velocity at axial location = "+ axial_location,'FontSize',15);
        colormap("turbo");
        xlabel('x','FontSize',15);
        ylabel('y','FontSize',15);
        colorbar;
        break
    end 

    
end 

