function [Xi_k] = function_Xi_0(time,Const,Config)

N_noeuds = Config.N_noeuds; 
Const.it_troncon = 1;

[DX,X_grille]=cheb(N_noeuds-1,Config.Li(1));  
X_k{1} = X_grille;

%Calculation of Xin X_1 of beam 1 
for it_x = 1:N_noeuds
    Xi(:,it_x) = [Config.K_0{1}(X_grille(it_x));0;0;1;Config.D_prime{1}(X_grille(it_x));Config.D{1}(X_grille(it_x))*Config.K_0{1}(X_grille(it_x))];
end
Xi_k{1} = Xi;

for it_tubes = 2 : Config.nb_tubes 
    % Calculation of X_k in X1
    Config.it_tubes = it_tubes;
    [T X_k{it_tubes}] = ode45(@(X_1,X_k) face_2_face_0(X_1,X_k,Const,Config),flip(X_grille),0); % Solve ODE
    
    % Calculation of Xi_k
    
    for it_x = 1:N_noeuds
        theta = Config.theta_0{it_tubes}(X_k{it_tubes}(it_x));
        Rx = [1 0 0;0 cos(theta) -sin(theta);0 sin(theta) cos(theta)];
        Rz = [-1, 0, 0;0, -1, 0;0, 0, 1];
        r1_k = [(-(Config.D{1}(X_grille(it_x))*eye(3)+Config.D{it_tubes}(X_k{it_tubes}(it_x))*Rx)*[0;1;0])];
        g1_k = [Rx*Rz,r1_k; 0 0 0 1];
        gk_1 = inv(g1_k);
        Rk_1 = gk_1(1:3,1:3);
        rk_1 = gk_1(1:3,4);
        
        [Phi,Phi_prime] = Base_Phi_prime((X_k{it_tubes}(it_x)),time,Const,Config);

        nb_modes_theta = Const.dim_base_q_theta/(Config.nb_tubes-1);
        
        mu_k = [1 0 0 0 0 Config.D{it_tubes}(X_k{it_tubes}(it_x))]';
        theta_prime_0_k = Config.theta_0_prime{it_tubes}(X_k{it_tubes}(it_x));
        theta_0_k = Config.theta_0{it_tubes}(X_k{it_tubes}(it_x));
        
        nu_k = [theta_prime_0_k,0,0,0,cos(theta_0_k)*Config.D_prime{1}(X_grille(it_x))-Config.D_prime{it_tubes}(X_k{it_tubes}(it_x)),sin(theta_0_k)*Config.D_prime{1}(X_grille(it_x))-Config.D{it_tubes}(X_k{it_tubes}(it_x))*theta_prime_0_k]';
        
        Xi_k_1 =  nu_k;
        
        Xi_k{it_tubes}(:,it_x) = Ad_g_(Rk_1,rk_1)*Xi(:,it_x) + Xi_k_1; 
    end
end