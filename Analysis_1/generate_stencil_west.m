
% This code generates a 2-D FVM stencil for western nodes.
% The faces of the cells do not have to be aligned to a 
% cartesian grid on any side, i.e. a cell can be any convex quadrilateral

% The code is divided in two parts. The first part generates the stencil
% and prints it to a file given in variable 'target file' after the 
% expression 'Start_stecil' (case sensivtive). The second part
% replaces the scalar expressions with matrix expression and allows to
% calculate the stencil for all inner nodes at once. This accelerates the
% calculation of the system Matrix a lot. 

% Either way all distances and surfaces (areas) have to be provided by you.

% Written by Thomas Runte


clear all 
clc

% Print results to file 
target_file = 'build_west.m';

fclose(fopen(target_file, 'w'));

%% First part
% Generate stencil with variables names introduced in Camilo F. Silva's
% Course "Numerical Thermo Fluid Dynamics"

% Initialize symbolic variables for distances,
% areas, lambdas, Temperatures

% Around eta_s
syms dy_S_Se dy_Se_e dy_e_P dy_P_S real
syms dx_S_Se dx_Se_e dx_e_P dx_P_S real

% Around eta_e
syms dy_s_sE dy_sE_nE dy_nE_n dy_n_s real
syms dx_s_sE dx_sE_nE dx_nE_n dx_n_s real

% Around eta_n
syms dy_P_e dy_e_Ne dy_Ne_N dy_N_P  real
syms dx_P_e dx_e_Ne dx_Ne_N dx_N_P  real 

% Around P
syms dy_s_se dy_se_ne dy_ne_n dy_n_s real
syms dx_s_se dx_se_ne dx_ne_n dx_n_s real

% Areas
syms S_eta S_eta_s S_e S_eta_n real

% Temperatures
syms  T_N T_P T_S T_NE T_E T_SE real

% inner Temperatures
syms T_eta T_eta_s T_eta_S T_eta_n T_eta_N T_se T_ne T_s T_n real

% Define inner Temperatures as interpolation of outer Temperatures
T_se=(T_SE+T_S+T_P+T_E)/4; %
T_ne=(T_NE+T_N+T_P+T_E)/4; %

T_e = (T_P + T_E)/2;
T_s = (T_P + T_S)/2;
T_n = (T_P + T_N)/2;

T_eta = (T_P + T_e)/2;
T_eta_s = (T_s + T_se)/2;
T_eta_n = (T_n + T_ne)/2;


% Gradients (Greens theorem)
dTdx_eta_s=  (dy_P_S*T_s + dy_S_Se*T_eta_s + dy_Se_e*T_se + dy_e_P*T_eta) /S_eta_s; %
dTdy_eta_s=  -(dx_P_S*T_s + dx_S_Se*T_eta_s + dx_Se_e*T_se + dx_e_P*T_eta) /S_eta_s; %

dTdx_e=  (dy_s_sE*T_se + dy_sE_nE*T_E + dy_nE_n*T_ne + dy_n_s*T_P) /S_e; %
dTdy_e=  -(dx_s_sE*T_se + dx_sE_nE*T_E + dx_nE_n*T_ne + dx_n_s*T_P) /S_e; %

dTdx_eta_n=  (dy_P_e*T_eta + dy_e_Ne*T_ne + dy_Ne_N*T_eta_n + dy_N_P*T_n) /S_eta_n; %
dTdy_eta_n=  -(dx_P_e*T_eta + dx_e_Ne*T_ne + dx_Ne_N*T_eta_n + dx_N_P*T_n) /S_eta_n; %


% Build hole stecil acounting for quadratic lambda like in Helmholtz 

 DDT= (dy_s_se*dTdx_eta_s - dx_s_se*dTdy_eta_s...
     + dy_se_ne*dTdx_e - dx_se_ne*dTdy_e...
     + dy_ne_n*dTdx_eta_n - dx_ne_n*dTdy_eta_n)/S_eta;
 
% Make Temperatur vector 
T=[T_E; T_S; T_N; T_NE; T_SE; T_P];


% NOTE: The Jacobian function is a "smart" way of factorizing the
% expresssion DDT with respect the temperature vector

stecil=jacobian(DDT,T);

% Find position in file
fileID2 = fopen(target_file, 'r+');



        fprintf(fileID2,'\n\n');
        fprintf(fileID2,'%% Nomenclature:\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%     N(i-1,j) -  Ne     NE(i-1,j+1)\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%                 |                 |\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%        n ------ ne - - - nE\n');
        fprintf(fileID2,'%%                 |                 |\n');
        fprintf(fileID2,'%%       |         |        |        |       |\n');
        fprintf(fileID2,'%%                 |                 |\n');
        fprintf(fileID2,'%%    P (i,j) - - e - -  E (i,j+1)\n');
        fprintf(fileID2,'%%                 |                 |\n');
        fprintf(fileID2,'%%       |         |        |        |       |\n');
        fprintf(fileID2,'%%                 |                 |\n');
        fprintf(fileID2,'%%      s ------ se - - - sE\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%                 |                 |\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%%   S(i+1,j)  - Se      SE(i+1,j+1)\n');
        fprintf(fileID2,'%%\n');
        fprintf(fileID2,'%% Indexing of stecil: \n\n');
        fprintf(fileID2,'%%   D_1 - D2\n');
        fprintf(fileID2,'%%     |     |     | \n');
        fprintf(fileID2,'%%   D_0 - D3\n');
        fprintf(fileID2,'%%     |     |     | \n');
        fprintf(fileID2,'%%    D1 - D4\n\n');

        
        %fprintf(fileID2,'lambda=boundary.lambda; \n\n');
      
        fprintf(fileID2,'%% Stecil \n\n');
        fprintf(fileID2,'%% East \n');
        fprintf(fileID2,'D3=%s; \n\n',char(stecil(1)));

        fprintf(fileID2,'%% South \n');
        fprintf(fileID2,'D1=%s; \n\n',char(stecil(2)));

        fprintf(fileID2,'%% North \n');
        fprintf(fileID2,'D_1=%s; \n\n',char(stecil(3)));

        fprintf(fileID2,'%% NE \n');
        fprintf(fileID2,'D2=%s; \n\n',char(stecil(4)));

        fprintf(fileID2,'%% SE \n');
        fprintf(fileID2,'D4=%s; \n\n',char(stecil(5)));

        fprintf(fileID2,'%% P \n');
        fprintf(fileID2,'D0=%s; \n\n',char(stecil(6)));

fclose(fileID2);
