% This function's objective is to apply the calculations through all the
% nozzle length. 
% The inputs are:
%           1. A_t: Throat area of the nozzle being considered
%           2. d_L: Axial location (in per one, so from 0= throat to 1=
%           exhaust)
%           3. A_e: Exhaust area of the nozzle being considered
% The outputs are:
%           1. Area_dif:   Area of cross section for each d_L
%           2. P_dif: pressure gradient insde coolant tube for each
%           d_L;


function [Area_dif, P_dif] = Complete_nozzle(D_t, d_L, D_e)

%% Initialization

pi = 3.1415;
A_t = pi*(D_t^2)/4;         % Throat Area 
A_e = pi*(D_e^2)/4;         % Exhaust Area
LTH = 1;                    %m
x_LTH = LTH*d_L;            %Point of interest

%% Geometry model: 
x = linspace(0,LTH);
y = ((D_e-D_t)/(2*log(2)))*(log(x+1)) +(D_t/2);
y_sym = -y;

y_dif = ((D_e-D_t)/(2*log(2)))*(log(x_LTH+1)) +(D_t/2);

% figure(1);
% line([x_LTH x_LTH],[-y_dif y_dif],'color','r','LineWidth',2);
% hold on;
% h(1)=plot(x,y,'-k','LineWidth',2,...
%     'MarkerSize',3,...
%     'MarkerEdgeColor','g');
% hold on;
% h(2)=plot(x,y_sym,'-k','LineWidth',2,...
%     'MarkerSize',3,...
%     'MarkerEdgeColor','g');
% grid on;
% title('Nozzle Geometry and analized axial location');
% hold off;

Area_dif =  pi*(y_dif^2);

P_dif_total = 100;          %from Throat (P=1) to Exhaust (P=100) f.e.
P_dif = P_dif_total*d_L;    %from Throat to point of interest.

