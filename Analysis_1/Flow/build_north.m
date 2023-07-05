

% Nomenclature:
%
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
% Indexing of stecil: 

%    D_3 -  D0 - D3
%     |     |     | 
%    D_2 -  D1 - D4

% Stecil 

% P 
D0=((dx_e_se*(dx_P_E/2 + (3*dx_s_P)/4 + dx_sE_s/4))/S_eta_e + (dx_sw_w*(dx_W_P/2 + (3*dx_P_s)/4 + dx_s_sW/4))/S_eta_w + (dy_e_se*(dy_P_E/2 + (3*dy_s_P)/4 + dy_sE_s/4))/S_eta_e + (dy_sw_w*(dy_W_P/2 + (3*dy_P_s)/4 + dy_s_sW/4))/S_eta_w + (dx_se_sw*(dx_w_e + dx_e_Se/4 + dx_Sw_w/4))/S_s + (dy_se_sw*(dy_w_e + dy_e_Se/4 + dy_Sw_w/4))/S_s)/S_eta; 

% East 
D3=((dx_e_se*(dx_P_E/2 + (3*dx_E_sE)/4 + dx_sE_s/4))/S_eta_e + (dy_e_se*(dy_P_E/2 + (3*dy_E_sE)/4 + dy_sE_s/4))/S_eta_e + (dx_e_Se*dx_se_sw)/(4*S_s) + (dy_e_Se*dy_se_sw)/(4*S_s))/S_eta; 

% SE 
D4=((dx_e_se*(dx_E_sE/4 + dx_sE_s/4))/S_eta_e + (dy_e_se*(dy_E_sE/4 + dy_sE_s/4))/S_eta_e + (dx_e_Se*dx_se_sw)/(4*S_s) + (dy_e_Se*dy_se_sw)/(4*S_s))/S_eta; 

% South 
D1=((dx_e_se*(dx_s_P/4 + dx_sE_s/4))/S_eta_e + (dx_sw_w*(dx_P_s/4 + dx_s_sW/4))/S_eta_w + (dy_e_se*(dy_s_P/4 + dy_sE_s/4))/S_eta_e + (dy_sw_w*(dy_P_s/4 + dy_s_sW/4))/S_eta_w + (dx_se_sw*(dx_e_Se/4 + dx_Sw_w/4 + dx_Se_Sw))/S_s + (dy_se_sw*(dy_e_Se/4 + dy_Sw_w/4 + dy_Se_Sw))/S_s)/S_eta; 

% SW 
D_2=((dx_sw_w*(dx_sW_W/4 + dx_s_sW/4))/S_eta_w + (dy_sw_w*(dy_sW_W/4 + dy_s_sW/4))/S_eta_w + (dx_Sw_w*dx_se_sw)/(4*S_s) + (dy_Sw_w*dy_se_sw)/(4*S_s))/S_eta; 

% West 
D_3=((dx_sw_w*(dx_W_P/2 + (3*dx_sW_W)/4 + dx_s_sW/4))/S_eta_w + (dy_sw_w*(dy_W_P/2 + (3*dy_sW_W)/4 + dy_s_sW/4))/S_eta_w + (dx_Sw_w*dx_se_sw)/(4*S_s) + (dy_Sw_w*dy_se_sw)/(4*S_s))/S_eta; 

