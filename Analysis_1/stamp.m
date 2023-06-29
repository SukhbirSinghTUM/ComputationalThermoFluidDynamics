function [stencil, B] = stamp(i, j, X, Y, B, TD)
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
%    verbose   verbositiy level
%
%  output:
%    stencil     linear equation for node (i,j)
%    b         new right-hand side value for node (i,j)


% Init

n = size(X, 1);
m = size(X, 2);
stencil = zeros(1, n*m);
index=@(ii, jj) ii + (jj-1)*n;


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

        build_inner

        %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
              
        % P
        stencil(index(i, j)) = D0;
        % East
        stencil(index(i, j+1)) = D3;
        % West
        stencil(index(i, j-1)) = D_3;
        % South
        stencil(index(i+1, j)) = D1;
        % North
        stencil(index(i-1, j)) = D_1;
        % NW
        stencil(index(i-1, j-1)) = D_4;
        % NE
        stencil(index(i-1, j+1)) = D2;
        % SW
        stencil(index(i+1, j-1)) = D_2;
        % SE
        stencil(index(i+1, j+1)) = D4;
        
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

        build_south

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
        B(index(i,j)) = -TD.flux*norm([dy_w_e dx_w_e])/S_eta;

    case 'North'

        % Principal node coordinates
        y_E = Y(i,j+1);      x_E = X(i,j+1);
        y_SE = Y(i+1,j+1);   x_SE = X(i+1,j+1);
        y_S = Y(i+1,j);      x_S = X(i+1,j);
        y_SW = Y(i+1,j-1);   x_SW = X(i+1,j-1);
        y_W = Y(i,j-1);      x_W = X(i,j-1);
        y_P = Y(i,j);        x_P = X(i,j);

        % Auxiliary node coordinates    
        y_sE = (y_E + y_SE)/2;  x_sE = (x_E + x_SE)/2;
        y_Se = (y_S + y_SE)/2;  x_Se = (x_S + x_SE)/2;
        y_Sw = (y_S + y_SW)/2;  x_Sw = (x_S + x_SW)/2;
        y_sW = (y_W + y_SW)/2;  x_sW = (x_W + x_SW)/2;

        y_w = (y_W + y_P)/2;  x_w = (x_W + x_P)/2;
        y_e = (y_E + y_P)/2;  x_e = (x_E + x_P)/2;
        y_s = (y_S + y_P)/2;  x_s = (x_S + x_P)/2;
       
        
        y_se = (y_s + y_sE)/2;  x_se = (x_s + x_sE)/2;
        y_sw = (y_s + y_sW)/2;  x_sw = (x_s + x_sW)/2;

        % Around s 
        dy_Se_Sw = y_Se - y_Sw;     dx_Se_Sw = x_Se - x_Sw;
        dy_e_Se = y_e - y_Se;       dx_e_Se = x_e - x_Se;
        dy_w_e = y_w - y_e;         dx_e_w = x_w - x_e;
        dy_Sw_w = y_Sw - y_w;       dx_Sw_w = x_Sw - x_w;
  
        % Around e 
        dy_sE_s = y_sE - y_s;       dx_sE_s = x_sE - x_s;
        dy_E_sE = y_E - y_sE;     dx_E_sE = x_E - x_sE;
        dy_P_E = y_P - y_E;       dx_P_E = x_P - x_E;
        dy_s_P = y_s - y_P;         dx_s_P = x_s - x_P;
        
        % Around w        
        dy_s_sW = y_s - y_sW;       dx_s_sW = x_s - x_sW;
        dy_P_s = y_P - y_s;         dx_P_s = x_P - x_s;
        dy_W_P = y_W - y_P;       dx_W_P = x_W - x_P;
        dy_sW_W = y_sW - y_W;     dx_sW_W = x_sW - x_W;

        % Around P 
        dy_se_sw = y_se - y_sw;     dx_se_sw = x_se - x_sw;
        dy_e_se = y_e - y_se;     dx_e_se = x_e - x_se;
        dy_w_e = y_w - y_e;     dx_w_e = x_w - x_e;
        dy_sw_w = y_sw - y_w;     dx_sw_w = x_sw - x_w;

        % Areas
        S_eta = abs((x_e*y_se - x_se*y_e)+(x_se*y_sw - x_sw*y_se)+(x_sw*y_w-x_w*y_sw)+(x_w*y_e - x_e*y_w))/2;
        S_s = abs((x_e*y_Se - x_Se*y_e)+(x_Se*y_Sw - x_Sw*y_Se)+(x_Sw*y_w - x_w*y_Sw)+(x_w*y_e - x_e*y_w))/2;
        S_eta_e = abs((x_E*y_sE - x_sE*y_E)+(x_sE*y_s - x_s*y_sE)+(x_s*y_P - x_P*y_s)+(x_P*y_E - x_E*y_P))/2;
        S_eta_w = abs((x_P*y_s - x_s*y_P)+(x_s*y_sW - x_sW*y_s)+(x_sW*y_W - x_W*y_sW)+(x_W*y_P - x_P*y_W))/2;

        %$$$$$$$$$$$$$$$$$$$$$$ Stencil $$$$$$$$$$$$$$$$$$$

        build_north

        %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

        % P
        stencil(index(i, j)) = D0 - (TD.htc*norm([dy_w_e dx_w_e]))/S_eta;
        % East
        stencil(index(i, j+1)) = D3;
        % West
        stencil(index(i, j-1)) = D_3;
        % South
        stencil(index(i+1, j)) = D1;
        % SW
        stencil(index(i+1, j-1)) = D_2;
        % SE
        stencil(index(i+1, j+1)) = D4;

        % Setting Boundary Conditions
        B(index(i,j)) = -TD.htc*TD.ambient*norm([dy_w_e dx_w_e])/S_eta;

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

        build_east

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

        % Principal node coordinates
        y_N  = Y(i-1,j);     x_N = X(i-1,j); 
        y_NE = Y(i-1,j+1);   x_NE = X(i-1,j+1);
        y_E = Y(i,j+1);      x_E = X(i,j+1);
        y_SE = Y(i+1,j+1);   x_SE = X(i+1,j+1);
        y_S = Y(i+1,j);      x_S = X(i+1,j);
        y_P = Y(i,j);        x_P = X(i,j);

        % Auxiliary node coordinates    
        y_Ne = (y_N + y_NE)/2;  x_Ne = (x_N + x_NE)/2;
        y_nE = (y_NE + y_E)/2;  x_nE = (x_NE + x_E)/2;
        y_sE = (y_E + y_SE)/2;  x_sE = (x_E + x_SE)/2;
        y_Se = (y_S + y_SE)/2;  x_Se = (x_S + x_SE)/2;

        y_e = (y_E + y_P)/2;  x_e = (x_E + x_P)/2;
        y_s = (y_S + y_P)/2;  x_s = (x_S + x_P)/2;
        y_n = (y_N + y_P)/2;  x_n = (x_N + x_P)/2;
       
        
        y_se = (y_s + y_sE)/2;  x_se = (x_s + x_sE)/2;
        y_ne = (y_n + y_nE)/2;  x_ne = (x_n + x_nE)/2;

        % Around s 
        dy_S_Se = y_Se - y_S;     dx_S_Se = x_Se - x_S;
        dy_Se_e = y_e - y_Se;       dx_Se_e = x_e - x_Se;
        dy_e_P = y_P - y_e;         dx_e_P = x_P - x_e;
        dy_P_S = y_S - y_P;       dx_P_S = x_S - x_P;
  
        % Around e
        dy_s_sE = y_sE - y_s;       dx_s_sE = x_sE - x_s;
        dy_sE_nE = y_nE - y_sE;     dx_sE_nE = x_nE - x_sE;
        dy_nE_n = y_n - y_nE;       dx_nE_n = x_n - x_nE;
        dy_n_s = y_s - y_n;         dx_n_s = x_s - x_n;

        % Around n 
        dy_P_e = y_e - y_P;         dx_P_e = x_e - x_P;
        dy_e_Ne = y_Ne - y_e;       dx_e_Ne = x_Ne - x_e;
        dy_Ne_N = y_N - y_Ne;     dx_Ne_N = x_N - x_Ne;
        dy_N_P = y_P - y_N;       dx_N_P = x_P - x_N;

        % Around P 
        dy_s_se = y_se - y_s;     dx_s_se = x_se - x_s;
        dy_se_ne = y_ne - y_se;     dx_se_ne = x_ne - x_se;
        dy_ne_n = y_n - y_ne;     dx_ne_n = x_n - x_ne;
        dy_n_s = y_s - y_n;     dx_n_s = x_s - x_n;

        % Areas
        S_eta = abs((x_ne*y_se - x_se*y_ne)+(x_se*y_s - x_s*y_se)+(x_s*y_n-x_n*y_s)+(x_n*y_ne - x_ne*y_n))/2;
        S_eta_s = abs((x_e*y_Se - x_Se*y_e)+(x_Se*y_S - x_S*y_Se)+(x_S*y_P - x_P*y_S)+(x_P*y_e - x_e*y_P))/2;
        S_e = abs((x_nE*y_sE - x_sE*y_nE)+(x_sE*y_s - x_s*y_sE)+(x_s*y_n - x_n*y_s)+(x_n*y_nE - x_nE*y_n))/2;
        S_eta_n = abs((x_Ne*y_e - x_e*y_Ne)+(x_e*y_P - x_P*y_e)+(x_P*y_N - x_N*y_P)+(x_N*y_Ne - x_Ne*y_N))/2;

        %$$$$$$$$$$$$$$$$$$$$$$ Stencil $$$$$$$$$$$$$$$$$$$

        build_west

        %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
              
        % P
        stencil(index(i, j)) = D0;
        % East
        stencil(index(i, j+1)) = D3;
        % South
        stencil(index(i+1, j)) = D1;
        % North
        stencil(index(i-1, j)) = D_1;
        % NE
        stencil(index(i-1, j+1)) = D2;
        % SE
        stencil(index(i+1, j+1)) = D4;

    case 'NorthWest'
        stencil(index(i, j)) = 1;
        stencil(index(i, j+1)) = -0.5;
        stencil(index(i+1, j)) = -0.5;

    case 'NorthEast'
        stencil(index(i, j)) = 1;
        stencil(index(i, j-1)) = -0.5;
        stencil(index(i+1, j)) = -0.5;

    case 'SouthWest'
        stencil(index(i, j)) = 1;
        stencil(index(i, j+1)) = -0.5;
        stencil(index(i-1, j)) = -0.5;

    case 'SouthEast'
        stencil(index(i, j)) = 1;
        stencil(index(i, j-1)) = -0.5;
        stencil(index(i-1, j)) = -0.5;

        
end



