function [stencil, B] = stamp_temp(i, j, X, Y, B, T_prev, u, Q_south, Tw_flow_solid2, Tw_west)
%stencil calculate the linear equation for node (i, j)

%  input:
%    i         node number in x direction
%    j         node number in y direction
%    X         x position of the nodes
%    Y         y position of the nodes
%    b         right-hand side value for node (i,j)
%    alpha     alpha
%    Tinf      Tinf for Robin BC
%    boundary  defines the boundary conditions
%
%  output:
%    stencil     linear equation for node (i,j)
%    b         new right-hand side value for node (i,j)


% Init
n = size(X, 1);
m = size(X, 2);
stencil = zeros(1, n*m);
index=@(ii, jj) ii + (jj-1)*n;

T_prev = reshape(T_prev,n*m,1);
u = reshape(u,n*m,1);

alpha = 10;                % Thermal diffusivity: alpha = 0.1m^2/s
delta_x = 0.1;              % axial step = 1cm


% Determine the node positon
if (index(i,j) <= n)&&(index(i,j)~=1)&&(index(i,j)~=n)
    nodePosition = 'West';
elseif (index(i,j)==1)
    nodePosition = 'NorthWest';
elseif (index(i,j)==n)
    nodePosition = 'SouthWest';
elseif (index(i,j)>n*(m-1)+1)&&(index(i,j)~=n*m)
    nodePosition = 'East';
elseif index(i,j)==(n*(m-1)+1)
    nodePosition = 'NorthEast';
elseif (index(i,j)==n*m)
    nodePosition = 'SouthEast';
elseif (rem(index(i,j),n)==0)&&(index(i,j)~=n)&&(index(i,j)~=n*m)
    nodePosition = 'South';
elseif (rem(index(i,j),n)==1)&&(index(i,j)~=1)&&(index(i,j)~=(n*(m-1)+1))
    nodePosition = 'North';
else 
    nodePosition = 'inner Node';
end



% Calculate the equation for the correct node position
switch nodePosition
    
    case 'inner Node'
        
% Nomenclature:
%
%    NW(i-1,j-1)   Nw -  N(i-1,j) -  Ne     NE(i-1,j+1)
%
%                 |                 |
%
%       nW - - - - nw ------ n ------ ne - - - nE
%                 |                 |
%       |         |        |        |       |
%                 |                 |
%   W(i, j-1) - - w - - P (i,j) - - e - -  E (i,j+1)
%                 |                 |
%       |         |        |        |       |
%                 |                 |
%      sW - - - - sw ------ s ------ se - - - sE
%
%                 |                 |
%
%   SW(i+1,j-1)   Sw  -  S(i+1,j)  - Se      SE(i+1,j+1)
%
% Indexing of stencil: 

%    D_4 - D_1 - D2
%     |     |     | 
%    D_3 - D_0 - D3
%     |     |     | 
%    D_2 -  D1 - D4
% Principal node coordinates
        y_NW = Y(i-1,j-1);   x_NW = X(i-1,j-1);
        y_N  = Y(i-1,j);     x_N = X(i-1,j); 
        y_NE = Y(i-1,j+1);   x_NE = X(i-1,j+1);
        y_E = Y(i,j+1);      x_E = X(i,j+1);
        y_SE = Y(i+1,j+1);   x_SE = X(i+1,j+1);
        y_S = Y(i+1,j);      x_S = X(i+1,j);
        y_SW = Y(i+1,j-1);   x_SW = X(i+1,j-1);
        y_W = Y(i,j-1);      x_W = X(i,j-1);
        y_P = Y(i,j);        x_P = X(i,j);

% Auxiliary node coordinates    
        y_Nw = (y_NW + y_N)/2;  x_Nw = (x_NW + x_N)/2; 
        y_Ne = (y_N + y_NE)/2;  x_Ne = (x_N + x_NE)/2;
        y_nE = (y_NE + y_E)/2;  x_nE = (x_NE + x_E)/2;
        y_sE = (y_E + y_SE)/2;  x_sE = (x_E + x_SE)/2;
        y_Se = (y_S + y_SE)/2;  x_Se = (x_S + x_SE)/2;
        y_Sw = (y_S + y_SW)/2;  x_Sw = (x_S + x_SW)/2;
        y_sW = (y_W + y_SW)/2;  x_sW = (x_W + x_SW)/2;
        y_nW = (y_W + y_NW)/2;  x_nW = (x_W + x_NW)/2;

        y_w = (y_W + y_P)/2;  x_w = (x_W + x_P)/2;
        y_e = (y_E + y_P)/2;  x_e = (x_E + x_P)/2;
        y_s = (y_S + y_P)/2;  x_s = (x_S + x_P)/2;
        y_n = (y_N + y_P)/2;  x_n = (x_N + x_P)/2;
       
        
        y_se = (y_s + y_sE)/2;  x_se = (x_s + x_sE)/2;
        y_sw = (y_s + y_sW)/2;  x_sw = (x_s + x_sW)/2;
        y_nw = (y_n + y_nW)/2;  x_nw = (x_n + x_nW)/2;
        y_ne = (y_n + y_nE)/2;  x_ne = (x_n + x_nE)/2;
  

        % Around s 
        dy_Sw_Se = y_Se - y_Sw;     dx_Sw_Se = x_Se - x_Sw;
        dy_Se_e = y_e - y_Se;       dx_Se_e = x_e - x_Se;
        dy_e_w = y_w - y_e;         dx_e_w = x_w - x_e;
        dy_w_Sw = y_Sw - y_w;       dx_w_Sw = x_Sw - x_w;
  
        % Around e
        dy_s_sE = y_sE - y_s;       dx_s_sE = x_sE - x_s;
        dy_sE_nE = y_nE - y_sE;     dx_sE_nE = x_nE - x_sE;
        dy_nE_n = y_n - y_nE;       dx_nE_n = x_n - x_nE;
        dy_n_s = y_s - y_n;         dx_n_s = x_s - x_n;

        % Around n 
        dy_w_e = y_e - y_w;         dx_w_e = x_e - x_w;
        dy_e_Ne = y_Ne - y_e;       dx_e_Ne = x_Ne - x_e;
        dy_Ne_Nw = y_Nw - y_Ne;     dx_Ne_Nw = x_Nw - x_Ne;
        dy_Nw_w = y_w - y_Nw;       dx_Nw_w = x_w - x_Nw;
        
        % Around w        
        dy_sW_s = y_s - y_sW;       dx_sW_s = x_s - x_sW;
        dy_s_n = y_n - y_s;         dx_s_n = x_n - x_s;
        dy_n_nW = y_nW - y_n;       dx_n_nW = x_nW - x_n;
        dy_nW_sW = y_sW - y_nW;     dx_nW_sW = x_sW - x_nW;

        % Around P
        dy_sw_se = y_se - y_sw;     dx_sw_se = x_se - x_sw;
        dy_se_ne = y_ne - y_se;     dx_se_ne = x_ne - x_se;
        dy_ne_nw = y_nw - y_ne;     dx_ne_nw = x_nw - x_ne;
        dy_nw_sw = y_sw - y_nw;     dx_nw_sw = x_sw - x_nw;

        % Areas
        S_P = abs((x_ne*y_se - x_se*y_ne)+(x_se*y_sw - x_sw*y_se)+(x_sw*y_nw-x_nw*y_sw)+(x_nw*y_ne - x_ne*y_nw))/2;
        S_s = abs((x_e*y_Se - x_Se*y_e)+(x_Se*y_Sw - x_Sw*y_Se)+(x_Sw*y_w - x_w*y_Sw)+(x_w*y_e - x_e*y_w))/2;
        S_e = abs((x_nE*y_sE - x_sE*y_nE)+(x_sE*y_s - x_s*y_sE)+(x_s*y_n - x_n*y_s)+(x_n*y_nE - x_nE*y_n))/2;
        S_n = abs((x_Ne*y_e - x_e*y_Ne)+(x_e*y_w - x_w*y_e)+(x_w*y_Nw - x_Nw*y_w)+(x_Nw*y_Ne - x_Ne*y_Nw))/2;
        S_w = abs((x_n*y_s - x_s*y_n)+(x_s*y_sW - x_sW*y_s)+(x_sW*y_nW - x_nW*y_sW)+(x_nW*y_n - x_n*y_nW))/2;
        
        %$$$$$$$$$$$$$$$$$$$$$$ Stencil $$$$$$$$$$$$$$$$$$$

        % East 
        D3=((dx_se_ne*(dx_nE_n/4 + dx_s_sE/4 + dx_sE_nE))/S_e + (dy_se_ne*(dy_nE_n/4 + dy_s_sE/4 + dy_sE_nE))/S_e + (dx_e_Ne*dx_ne_nw)/(4*S_n) + (dx_Se_e*dx_sw_se)/(4*S_s) + (dy_e_Ne*dy_ne_nw)/(4*S_n) + (dy_Se_e*dy_sw_se)/(4*S_s))/S_P; 
        
        % West 
        D_3=((dx_nw_sw*(dx_n_nW/4 + dx_sW_s/4 + dx_nW_sW))/S_w + (dy_nw_sw*(dy_n_nW/4 + dy_sW_s/4 + dy_nW_sW))/S_w + (dx_Nw_w*dx_ne_nw)/(4*S_n) + (dx_w_Sw*dx_sw_se)/(4*S_s) + (dy_Nw_w*dy_ne_nw)/(4*S_n) + (dy_w_Sw*dy_sw_se)/(4*S_s))/S_P; 
        
        % South 
        D1=((dx_sw_se*(dx_Se_e/4 + dx_w_Sw/4 + dx_Sw_Se))/S_s + (dy_sw_se*(dy_Se_e/4 + dy_w_Sw/4 + dy_Sw_Se))/S_s + (dx_s_sE*dx_se_ne)/(4*S_e) + (dx_sW_s*dx_nw_sw)/(4*S_w) + (dy_s_sE*dy_se_ne)/(4*S_e) + (dy_sW_s*dy_nw_sw)/(4*S_w))/S_P; 
        
        % North 
        D_1=((dx_ne_nw*(dx_e_Ne/4 + dx_Nw_w/4 + dx_Ne_Nw))/S_n + (dy_ne_nw*(dy_e_Ne/4 + dy_Nw_w/4 + dy_Ne_Nw))/S_n + (dx_nE_n*dx_se_ne)/(4*S_e) + (dx_n_nW*dx_nw_sw)/(4*S_w) + (dy_nE_n*dy_se_ne)/(4*S_e) + (dy_n_nW*dy_nw_sw)/(4*S_w))/S_P; 
        
        % NW 
        D_4=((dx_Nw_w*dx_ne_nw)/(4*S_n) + (dx_n_nW*dx_nw_sw)/(4*S_w) + (dy_Nw_w*dy_ne_nw)/(4*S_n) + (dy_n_nW*dy_nw_sw)/(4*S_w))/S_P; 
        
        % NE 
        D2=((dx_nE_n*dx_se_ne)/(4*S_e) + (dx_e_Ne*dx_ne_nw)/(4*S_n) + (dy_nE_n*dy_se_ne)/(4*S_e) + (dy_e_Ne*dy_ne_nw)/(4*S_n))/S_P; 
        
        % SW 
        D_2=((dx_w_Sw*dx_sw_se)/(4*S_s) + (dx_sW_s*dx_nw_sw)/(4*S_w) + (dy_w_Sw*dy_sw_se)/(4*S_s) + (dy_sW_s*dy_nw_sw)/(4*S_w))/S_P; 
        
        % SE 
        D4=((dx_s_sE*dx_se_ne)/(4*S_e) + (dx_Se_e*dx_sw_se)/(4*S_s) + (dy_s_sE*dy_se_ne)/(4*S_e) + (dy_Se_e*dy_sw_se)/(4*S_s))/S_P; 
        
        % P 
        D0=((dx_se_ne*(dx_n_s + dx_nE_n/4 + dx_s_sE/4))/S_e + (dx_ne_nw*(dx_w_e + dx_e_Ne/4 + dx_Nw_w/4))/S_n + (dx_sw_se*(dx_e_w + dx_Se_e/4 + dx_w_Sw/4))/S_s + (dx_nw_sw*(dx_s_n + dx_n_nW/4 + dx_sW_s/4))/S_w + (dy_se_ne*(dy_n_s + dy_nE_n/4 + dy_s_sE/4))/S_e + (dy_ne_nw*(dy_w_e + dy_e_Ne/4 + dy_Nw_w/4))/S_n + (dy_sw_se*(dy_e_w + dy_Se_e/4 + dy_w_Sw/4))/S_s + (dy_nw_sw*(dy_s_n + dy_n_nW/4 + dy_sW_s/4))/S_w)/S_P;

        %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
              
        % P
        stencil(index(i, j)) = (alpha*delta_x/u(index(i,j)))*D0 - 1;
        % East
        stencil(index(i, j+1)) = (alpha*delta_x/u(index(i,j)))*D3;
        % West
        stencil(index(i, j-1)) = (alpha*delta_x/u(index(i,j)))*D_3;
        % South
        stencil(index(i+1, j)) = (alpha*delta_x/u(index(i,j)))*D1;
        % North
        stencil(index(i-1, j)) = (alpha*delta_x/u(index(i,j)))*D_1;
        % NW
        stencil(index(i-1, j-1)) = (alpha*delta_x/u(index(i,j)))*D_4;
        % NE
        stencil(index(i-1, j+1)) = (alpha*delta_x/u(index(i,j)))*D2;
        % SW
        stencil(index(i+1, j-1)) = (alpha*delta_x/u(index(i,j)))*D_2;
        % SE
        stencil(index(i+1, j+1)) = (alpha*delta_x/u(index(i,j)))*D4;

        % Boundary Condition
        B(index(i,j)) = -T_prev(index(i,j));


    case 'North'
        stencil(index(i,j)) = 1;
        B(index(i,j)) = Tw_flow_solid2(j);

    case 'South'
        % Principal node coordinates
        y_NW = Y(i-1,j-1);   x_NW = X(i-1,j-1);
        y_N  = Y(i-1,j);     x_N = X(i-1,j); 
        y_NE = Y(i-1,j+1);   x_NE = X(i-1,j+1);
        y_E = Y(i,j+1);      x_E = X(i,j+1);
        y_W = Y(i,j-1);      x_W = X(i,j-1);
        y_P = Y(i,j);        x_P = X(i,j);

        % Auxiliary node coordinates    
        y_Nw = (y_NW + y_N)/2;  x_Nw = (x_NW + x_N)/2; 
        y_Ne = (y_N + y_NE)/2;  x_Ne = (x_N + x_NE)/2;
        y_nE = (y_NE + y_E)/2;  x_nE = (x_NE + x_E)/2;
        y_nW = (y_W + y_NW)/2;  x_nW = (x_W + x_NW)/2;

        y_w = (y_W + y_P)/2;  x_w = (x_W + x_P)/2;
        y_e = (y_E + y_P)/2;  x_e = (x_E + x_P)/2;
        y_n = (y_N + y_P)/2;  x_n = (x_N + x_P)/2;
       
        y_nw = (y_n + y_nW)/2;  x_nw = (x_n + x_nW)/2;
        y_ne = (y_n + y_nE)/2;  x_ne = (x_n + x_nE)/2;

        % Distances around P
        dy_w_e = y_e - y_w;         dx_w_e = x_e - x_w;
        dy_e_ne = y_ne - y_e;       dx_e_ne = x_ne - x_e;
        dy_ne_nw = y_nw - y_ne;     dx_ne_nw = x_nw - x_ne;
        dy_nw_w = y_w - y_nw;       dx_nw_w = x_w - x_nw;

        % Around e
        dy_P_E = y_E - y_P;         dx_P_E = x_E - x_P;
        dy_E_nE = y_nE - y_E;       dx_E_nE = x_nE - x_E;
        dy_nE_n = y_n - y_nE;       dx_nE_n = x_n - x_nE;
        dy_n_P = y_P - y_n;         dx_n_P = x_P - x_n;

        % Around w        
        dy_W_P = y_P - y_W;         dx_W_P = x_P - x_W;
        dy_P_n = y_n - y_P;         dx_P_n = x_n - x_P;
        dy_n_nW = y_nW - y_n;       dx_n_nW = x_nW - x_n;
        dy_nW_W = y_W - y_nW;       dx_nW_W = x_W - x_nW;

        % Around n 
        dy_w_e = y_e - y_w;         dx_w_e = x_e - x_w;
        dy_e_Ne = y_Ne - y_e;       dx_e_Ne = x_Ne - x_e;
        dy_Ne_Nw = y_Nw - y_Ne;     dx_Ne_Nw = x_Nw - x_Ne;
        dy_Nw_w = y_w - y_Nw;       dx_Nw_w = x_w - x_Nw;

        % Areas
        S_eta = abs((x_ne*y_e - x_e*y_ne)+(x_e*y_w - x_w*y_e)+(x_w*y_nw-x_nw*y_w)+(x_nw*y_ne - x_ne*y_nw))/2;
        S_eta_e = abs((x_nE*y_E - x_E*y_nE)+(x_E*y_P - x_P*y_E)+(x_P*y_n - x_n*y_P)+(x_n*y_nE - x_nE*y_n))/2;
        S_n = abs((x_Ne*y_e - x_e*y_Ne)+(x_e*y_w - x_w*y_e)+(x_w*y_Nw - x_Nw*y_w)+(x_Nw*y_Ne - x_Ne*y_Nw))/2;
        S_eta_w = abs((x_n*y_P - x_P*y_n)+(x_P*y_W - x_W*y_P)+(x_W*y_nW - x_nW*y_W)+(x_nW*y_n - x_n*y_nW))/2;

        %$$$$$$$$$$$$$$$$$$$$$$ Stencil $$$$$$$$$$$$$$$$$$$

        % P 
        D0=((dx_e_ne*(dx_P_E/2 + (3*dx_n_P)/4 + dx_nE_n/4))/S_eta_e + (dx_nw_w*(dx_W_P/2 + (3*dx_P_n)/4 + dx_n_nW/4))/S_eta_w + (dy_e_ne*(dy_P_E/2 + (3*dy_n_P)/4 + dy_nE_n/4))/S_eta_e + (dy_nw_w*(dy_W_P/2 + (3*dy_P_n)/4 + dy_n_nW/4))/S_eta_w + (dx_ne_nw*(dx_w_e + dx_e_Ne/4 + dx_Nw_w/4))/S_n + (dy_ne_nw*(dy_w_e + dy_e_Ne/4 + dy_Nw_w/4))/S_n)/S_eta; 
        
        % East 
        D3=((dx_e_ne*(dx_P_E/2 + (3*dx_E_nE)/4 + dx_nE_n/4))/S_eta_e + (dy_e_ne*(dy_P_E/2 + (3*dy_E_nE)/4 + dy_nE_n/4))/S_eta_e + (dx_e_Ne*dx_ne_nw)/(4*S_n) + (dy_e_Ne*dy_ne_nw)/(4*S_n))/S_eta; 
        
        % NE 
        D2=((dx_e_ne*(dx_E_nE/4 + dx_nE_n/4))/S_eta_e + (dy_e_ne*(dy_E_nE/4 + dy_nE_n/4))/S_eta_e + (dx_e_Ne*dx_ne_nw)/(4*S_n) + (dy_e_Ne*dy_ne_nw)/(4*S_n))/S_eta; 
        
        % North 
        D_1=((dx_e_ne*(dx_n_P/4 + dx_nE_n/4))/S_eta_e + (dx_nw_w*(dx_P_n/4 + dx_n_nW/4))/S_eta_w + (dy_e_ne*(dy_n_P/4 + dy_nE_n/4))/S_eta_e + (dy_nw_w*(dy_P_n/4 + dy_n_nW/4))/S_eta_w + (dx_ne_nw*(dx_e_Ne/4 + dx_Nw_w/4 + dx_Ne_Nw))/S_n + (dy_ne_nw*(dy_e_Ne/4 + dy_Nw_w/4 + dy_Ne_Nw))/S_n)/S_eta; 
        
        % NW 
        D_4=((dx_nw_w*(dx_nW_W/4 + dx_n_nW/4))/S_eta_w + (dy_nw_w*(dy_nW_W/4 + dy_n_nW/4))/S_eta_w + (dx_Nw_w*dx_ne_nw)/(4*S_n) + (dy_Nw_w*dy_ne_nw)/(4*S_n))/S_eta; 
        
        % West 
        D_3=((dx_nw_w*(dx_W_P/2 + (3*dx_nW_W)/4 + dx_n_nW/4))/S_eta_w + (dy_nw_w*(dy_W_P/2 + (3*dy_nW_W)/4 + dy_n_nW/4))/S_eta_w + (dx_Nw_w*dx_ne_nw)/(4*S_n) + (dy_Nw_w*dy_ne_nw)/(4*S_n))/S_eta; 

        %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
              
        % P
        stencil(index(i, j)) = D0;
        % East
        stencil(index(i, j+1)) = D3;
        % West
        stencil(index(i, j-1)) = D_3;
        % North
        stencil(index(i-1, j)) = D_1;
        % NW
        stencil(index(i-1, j-1)) = D_4;
        % NE
        stencil(index(i-1, j+1)) = D2;

        % Boundary Condition
        B(index(i,j)) = -Q_south(1,j)*norm([dy_w_e dx_w_e])/S_eta;
        

    case 'East'
        % Principal node coordinates
        y_NW = Y(i-1,j-1);   x_NW = X(i-1,j-1);
        y_N  = Y(i-1,j);     x_N = X(i-1,j); 
        y_S = Y(i+1,j);      x_S = X(i+1,j);
        y_SW = Y(i+1,j-1);   x_SW = X(i+1,j-1);
        y_W = Y(i,j-1);      x_W = X(i,j-1);
        y_P = Y(i,j);        x_P = X(i,j);

        % Auxiliary node coordinates    
        y_Nw = (y_NW + y_N)/2;  x_Nw = (x_NW + x_N)/2; 
        y_Sw = (y_S + y_SW)/2;  x_Sw = (x_S + x_SW)/2;
        y_sW = (y_W + y_SW)/2;  x_sW = (x_W + x_SW)/2;
        y_nW = (y_W + y_NW)/2;  x_nW = (x_W + x_NW)/2;

        y_w = (y_W + y_P)/2;  x_w = (x_W + x_P)/2;
        y_s = (y_S + y_P)/2;  x_s = (x_S + x_P)/2;
        y_n = (y_N + y_P)/2;  x_n = (x_N + x_P)/2;
      
        y_sw = (y_s + y_sW)/2;  x_sw = (x_s + x_sW)/2;
        y_nw = (y_n + y_nW)/2;  x_nw = (x_n + x_nW)/2;

        % Around s 
        dy_Sw_S = y_S - y_Sw;       dx_Sw_S = x_S - x_Sw;
        dy_S_P = y_P - y_S;         dx_S_P = x_P - x_S;
        dy_P_w = y_w - y_P;         dx_P_w = x_w - x_P;
        dy_w_Sw = y_Sw - y_w;       dx_w_Sw = x_Sw - x_w;

        % Around n 
        dy_w_P = y_P - y_w;         dx_w_P = x_P - x_w;
        dy_P_N = y_N - y_P;         dx_P_N = x_N - x_P;
        dy_N_Nw = y_Nw - y_N;       dx_N_Nw = x_Nw - x_N;
        dy_Nw_w = y_w - y_Nw;       dx_Nw_w = x_w - x_Nw;
        
        % Around w d    
        dy_sW_s = y_s - y_sW;       dx_sW_s = x_s - x_sW;
        dy_s_n = y_n - y_s;         dx_s_n = x_n - x_s;
        dy_n_nW = y_nW - y_n;       dx_n_nW = x_nW - x_n;
        dy_nW_sW = y_sW - y_nW;     dx_nW_sW = x_sW - x_nW;

        % Around P 
        dy_sw_s = y_s - y_sw;       dx_sw_s = x_s - x_sw;
        dy_s_n = y_n - y_s;         dx_s_n = x_n - x_s;
        dy_n_nw = y_nw - y_n;       dx_n_nw = x_nw - x_n;
        dy_nw_sw = y_sw - y_nw;     dx_nw_sw = x_sw - x_nw;

        % Areas
        S_eta = abs((x_n*y_s - x_s*y_n)+(x_s*y_sw - x_sw*y_s)+(x_sw*y_nw-x_nw*y_sw)+(x_nw*y_n - x_n*y_nw))/2;
        S_s_eta = abs((x_P*y_S - x_S*y_P)+(x_S*y_Sw - x_Sw*y_S)+(x_Sw*y_w - x_w*y_Sw)+(x_w*y_P - x_P*y_w))/2;
        S_n_eta = abs((x_N*y_P - x_P*y_N)+(x_P*y_w - x_w*y_P)+(x_w*y_Nw - x_Nw*y_w)+(x_Nw*y_N - x_N*y_Nw))/2;
        S_w = abs((x_n*y_s - x_s*y_n)+(x_s*y_sW - x_sW*y_s)+(x_sW*y_nW - x_nW*y_sW)+(x_nW*y_n - x_n*y_nW))/2;

        %$$$$$$$$$$$$$$$$$$$$$$ Stencil $$$$$$$$$$$$$$$$$$$

        % P 
        D0=((dx_n_nw*(dx_P_N/2 + (3*dx_w_P)/4 + dx_Nw_w/4))/S_n_eta + (dx_sw_s*(dx_S_P/2 + (3*dx_P_w)/4 + dx_w_Sw/4))/S_s_eta + (dy_n_nw*(dy_P_N/2 + (3*dy_w_P)/4 + dy_Nw_w/4))/S_n_eta + (dy_sw_s*(dy_S_P/2 + (3*dy_P_w)/4 + dy_w_Sw/4))/S_s_eta + (dx_nw_sw*(dx_s_n + dx_n_nW/4 + dx_sW_s/4))/S_w + (dy_nw_sw*(dy_s_n + dy_n_nW/4 + dy_sW_s/4))/S_w)/S_eta; 
        
        % North 
        D_1=((dx_n_nw*(dx_P_N/2 + (3*dx_N_Nw)/4 + dx_Nw_w/4))/S_n_eta + (dy_n_nw*(dy_P_N/2 + (3*dy_N_Nw)/4 + dy_Nw_w/4))/S_n_eta + (dx_n_nW*dx_nw_sw)/(4*S_w) + (dy_n_nW*dy_nw_sw)/(4*S_w))/S_eta; 
        
        % NW 
        D_4=((dx_n_nw*(dx_N_Nw/4 + dx_Nw_w/4))/S_n_eta + (dy_n_nw*(dy_N_Nw/4 + dy_Nw_w/4))/S_n_eta + (dx_n_nW*dx_nw_sw)/(4*S_w) + (dy_n_nW*dy_nw_sw)/(4*S_w))/S_eta; 
        
        % West 
        D_3=((dx_n_nw*(dx_w_P/4 + dx_Nw_w/4))/S_n_eta + (dx_sw_s*(dx_P_w/4 + dx_w_Sw/4))/S_s_eta + (dy_n_nw*(dy_w_P/4 + dy_Nw_w/4))/S_n_eta + (dy_sw_s*(dy_P_w/4 + dy_w_Sw/4))/S_s_eta + (dx_nw_sw*(dx_n_nW/4 + dx_sW_s/4 + dx_nW_sW))/S_w + (dy_nw_sw*(dy_n_nW/4 + dy_sW_s/4 + dy_nW_sW))/S_w)/S_eta; 
        
        % SW 
        D_2=((dx_sw_s*(dx_Sw_S/4 + dx_w_Sw/4))/S_s_eta + (dy_sw_s*(dy_Sw_S/4 + dy_w_Sw/4))/S_s_eta + (dx_sW_s*dx_nw_sw)/(4*S_w) + (dy_sW_s*dy_nw_sw)/(4*S_w))/S_eta; 
        
        % South 
        D1=((dx_sw_s*(dx_S_P/2 + (3*dx_Sw_S)/4 + dx_w_Sw/4))/S_s_eta + (dy_sw_s*(dy_S_P/2 + (3*dy_Sw_S)/4 + dy_w_Sw/4))/S_s_eta + (dx_sW_s*dx_nw_sw)/(4*S_w) + (dy_sW_s*dy_nw_sw)/(4*S_w))/S_eta; 

        %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

        % P
        stencil(index(i, j)) = D0;
        % West
        stencil(index(i, j-1)) = D_3;
        % South
        stencil(index(i+1, j)) = D1;
        % North
        stencil(index(i-1, j)) = D_1;
        % NW
        stencil(index(i-1, j-1)) = D_4;
        % SW
        stencil(index(i+1, j-1)) = D_2;

    case 'West'
        stencil(index(i,j)) = 1;
        B(index(i,j)) = Tw_west(i,1);

    case 'NorthWest'
        stencil(index(i, j)) = 1;
        stencil(index(i, j+1)) = -1;
        stencil(index(i+1, j)) = -1;
        stencil(index(i+1, j+1)) = 1;

    case 'NorthEast'
        stencil(index(i, j)) = 1;
        stencil(index(i, j-1)) = -1;
        stencil(index(i+1, j)) = -1;
        stencil(index(i+1, j-1)) = 1;

    case 'SouthWest'
        stencil(index(i, j)) = 1;
        stencil(index(i, j+1)) = -1;
        stencil(index(i-1, j)) = -1;
        stencil(index(i-1, j+1)) = 1;

    case 'SouthEast'
        stencil(index(i, j)) = 1;
        stencil(index(i, j-1)) = -1;
        stencil(index(i-1, j)) = -1;
        stencil(index(i-1, j-1)) = 1;
        
end



