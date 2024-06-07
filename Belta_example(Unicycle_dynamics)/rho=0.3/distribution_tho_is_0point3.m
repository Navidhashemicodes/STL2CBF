clear
clc
close all
% load('Result_rho_is_0point3.mat')
load('Result_rho_is_0point3_extendedto_40000.mat')

Param=param{iter};
Lam=lam_param{iter};
Bet=beta_param{iter};

NN=100000;
[JP, J_STL_h, J_STL_s] = expectations_Unicycle([0.6 1.4],[0.6 1.4],[pi/2-pi/10 , pi/2+pi/10],[0.99, 1.01], NN, M, P, Param, Lam, Bet, T);

