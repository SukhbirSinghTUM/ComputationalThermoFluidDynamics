

% Nomenclature:
%
%    NW(i-1,j-1)   Nw -  N(i-1,j)
%
%                 |              
%
%       nW - - - - nw ------ n   
%                 |              
%       |         |        |     
%                 |              
%   W(i, j-1) - - w - - P (i,j)  
%                 |              
%       |         |        |     
%                 |              
%      sW - - - - sw ------ s    
%
%                 |              
%
%   SW(i+1,j-1)   Sw  -  S(i+1,j)
%
% Indexing of stecil: 

%    D_4 - D_1 
%     |     |  
%    D_3 -  D0 
%     |     |  
%    D_2 -  D1 

% Stecil 

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

