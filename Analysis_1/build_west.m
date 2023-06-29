

% Nomenclature:
%
%     N(i-1,j) -  Ne     NE(i-1,j+1)
%
%                 |                 |
%
%        n ------ ne - - - nE
%                 |                 |
%       |         |        |        |       |
%                 |                 |
%    P (i,j) - - e - -  E (i,j+1)
%                 |                 |
%       |         |        |        |       |
%                 |                 |
%      s ------ se - - - sE
%
%                 |                 |
%
%   S(i+1,j)  - Se      SE(i+1,j+1)
%
% Indexing of stecil: 

%   D_1 - D2
%     |     |     | 
%   D_0 - D3
%     |     |     | 
%    D1 - D4

% Stecil 

% East 
D3=((dx_ne_n*(dx_P_e/4 + dx_Ne_N/8 + dx_e_Ne/4))/S_eta_n + (dx_s_se*(dx_e_P/4 + dx_S_Se/8 + dx_Se_e/4))/S_eta_s + (dy_ne_n*(dy_P_e/4 + dy_Ne_N/8 + dy_e_Ne/4))/S_eta_n + (dy_s_se*(dy_e_P/4 + dy_S_Se/8 + dy_Se_e/4))/S_eta_s + (dx_se_ne*(dx_nE_n/4 + dx_s_sE/4 + dx_sE_nE))/S_e + (dy_se_ne*(dy_nE_n/4 + dy_s_sE/4 + dy_sE_nE))/S_e)/S_eta; 

% South 
D1=((dx_s_se*(dx_P_S/2 + (3*dx_S_Se)/8 + dx_Se_e/4))/S_eta_s + (dy_s_se*(dy_P_S/2 + (3*dy_S_Se)/8 + dy_Se_e/4))/S_eta_s + (dx_s_sE*dx_se_ne)/(4*S_e) + (dy_s_sE*dy_se_ne)/(4*S_e))/S_eta; 

% North 
D_1=((dx_ne_n*(dx_N_P/2 + (3*dx_Ne_N)/8 + dx_e_Ne/4))/S_eta_n + (dy_ne_n*(dy_N_P/2 + (3*dy_Ne_N)/8 + dy_e_Ne/4))/S_eta_n + (dx_nE_n*dx_se_ne)/(4*S_e) + (dy_nE_n*dy_se_ne)/(4*S_e))/S_eta; 

% NE 
D2=((dx_ne_n*(dx_Ne_N/8 + dx_e_Ne/4))/S_eta_n + (dy_ne_n*(dy_Ne_N/8 + dy_e_Ne/4))/S_eta_n + (dx_nE_n*dx_se_ne)/(4*S_e) + (dy_nE_n*dy_se_ne)/(4*S_e))/S_eta; 

% SE 
D4=((dx_s_se*(dx_S_Se/8 + dx_Se_e/4))/S_eta_s + (dy_s_se*(dy_S_Se/8 + dy_Se_e/4))/S_eta_s + (dx_s_sE*dx_se_ne)/(4*S_e) + (dy_s_sE*dy_se_ne)/(4*S_e))/S_eta; 

% P 
D0=((dx_se_ne*(dx_n_s + dx_nE_n/4 + dx_s_sE/4))/S_e + (dy_se_ne*(dy_n_s + dy_nE_n/4 + dy_s_sE/4))/S_e + (dx_ne_n*(dx_N_P/2 + (3*dx_P_e)/4 + (3*dx_Ne_N)/8 + dx_e_Ne/4))/S_eta_n + (dx_s_se*(dx_P_S/2 + (3*dx_e_P)/4 + (3*dx_S_Se)/8 + dx_Se_e/4))/S_eta_s + (dy_ne_n*(dy_N_P/2 + (3*dy_P_e)/4 + (3*dy_Ne_N)/8 + dy_e_Ne/4))/S_eta_n + (dy_s_se*(dy_P_S/2 + (3*dy_e_P)/4 + (3*dy_S_Se)/8 + dy_Se_e/4))/S_eta_s)/S_eta; 

