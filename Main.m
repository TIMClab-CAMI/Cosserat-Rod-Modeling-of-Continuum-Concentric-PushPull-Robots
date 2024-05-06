clc
clear all
close all


% add path tool
addpath('./tool')
addpath('./spectral_method')

% numbers of nodes for spectral method
Config.N_noeuds = 40; 
Config.time  = 0;

% number of tube
Config.nb_tubes = 2;

% Lenght of the CCPPR
Config.Li = 150e-3;
Config.L = Config.Li;

if Config.nb_tubes == 3
    
    % definition of geometric variables
    Config.K_0{1} = @(X) 0*2*pi/(Config.Li);

    Config.D{1} = @(X) (1e-3)*(1.12);
    Config.D{2} = @(X) (1e-3)*(0.983);
    Config.D{3} = @(X) (1e-3)*(0.764);

    Config.D_prime{1} = @(X) (1e-3)*(0);
    Config.D_prime{2} = @(X) (1e-3)*(0);
    Config.D_prime{3} = @(X) (1e-3)*(0);

    Config.theta_0{1} = @(X) (0);
    Config.theta_0{2} = @(X) (pi/3);
    Config.theta_0{3} = @(X) (-pi/3);

    Config.theta_0_prime{1} = @(X) 0*(0);
    Config.theta_0_prime{2} = @(X) 0*(pi/3);
    Config.theta_0_prime{3} = @(X) 0*(-pi/3);
end

if Config.nb_tubes == 2
    
    % definition of geometric variables
    Config.K_0{1} = @(X) 2*pi/(Config.Li);

    Config.D{1} = @(X) (1e-3)*(3.024);
    Config.D{2} = @(X) (1e-3)*(1.134);
   
    Config.D_prime{1} = @(X) (1e-3)*(0);
    Config.D_prime{2} = @(X) (1e-3)*(0);

    Config.theta_0{1} = @(X) (0);
    Config.theta_0{2} = @(X) 0*(pi/3);

    Config.theta_0_prime{1} = @(X) 0*(0);
    Config.theta_0_prime{2} = @(X) 0*(pi/3);
end

% Definition of actuated variables reference tube
% K1, K2, K3, ...
% if value is 1 -> variable is activated
% if value is 0 -> variable not activated
Config.V_a = [1,1,1,0,0,0];

% Definition of the number of modes for each variable
Const.dim_base_k = [2,3,3,0,0,0];

% Size of q_epsilon (q_1)
Const.dim_base_q_e     = Config.V_a*Const.dim_base_k';
% Size of theta_k
Const.dim_base_q_theta = (Config.nb_tubes-1)*2;

% Size of q
Const.dim_base = Const.dim_base_q_e + Const.dim_base_q_theta;

% Gravity
Const.Gamma_g = 9.81;

% Physical parameters of the beam 
parametre_constant;

%% ------------------------------------------------------- %
%
% pseudo "temporal" integration

% Initialisation of q
Const.q         = zeros(Const.dim_base,1);

% Initial position and orientation of the base of the CCPPR
r_0 = [0;Config.D{1}(0);0];
q_0 = [1 0 0 0]';

Const.r_0 = r_0;
Const.q_0 = q_0; 

% Force values F- -> Fp (-->m)
Const.Fm_materielle = zeros(6,1); % Material
Const.Fm_spaciale = zeros(6,1); % Spacial

% Force values F+ -> Fp (+->p)
Const.Fp_materielle = zeros(6,1); % Material
Const.Fp_spaciale = 0*[0;0;0;-2*0.1327;0.3693;0];% Spacial

% Force values Fbar
Const.Fbar_materielle = zeros(6,1); % Material
Const.Fbar_spaciale = 0*[0;0;0;0;-Const.Gamma_g*Const.Aire*Const.rho;0]; % Spacial

% Initialisation of T and C
for it_tubes = 2 : Config.nb_tubes  
    Const.T{it_tubes} = 0;
    Const.C{it_tubes} = 0;
end

tic
[T,q] = Time_Integration_Newton_Static_Notch_Spectral(Config.time,Const,Config);
tsimul = toc
save('Test_Static_notch_gravite.mat');
