function [Q_a_T] = ISM_spectral_notch(q_0,r_0,q,a,b,time,Const,Config)
% Hookean value

if Config.nb_tubes == 3
    Hook{1} = @(X) (1e-4)*diag([2.73,1.78,1.96,1000,1000,1000]);%1.25e-4;%7e-4;%
    Hook{2} = @(X) (1e-4)*diag([2.78,1.78,1.81,1000,1000,1000]);%1.85e-4;%3e-4;% 
    Hook{3} = @(X) (1e-4)*diag([2.91,1.89,1.83,1000,1000,1000]);%1.85e-4;%3e-4;% 
end
% Hookean value
if Config.nb_tubes == 3
    Hook{1} = @(X) (1e-4)*diag([2.73,1.78,1.96,1000,1000,1000]);
    Hook{2} = @(X) (1e-4)*diag([2.78,1.78,1.81,1000,1000,1000]); 
    Hook{3} = @(X) (1e-4)*diag([2.91,1.89,1.83,1000,1000,1000]); 
end
if Config.nb_tubes == 2
    Hook{1} = @(X) diag([0.0039,0.0096*1.3,0.0006*1.75,9.2e3*1.75,3.4e3,3.4e3]);%1.25e-4;%7e-4;%
    Hook{2} = @(X) diag([0.0012,0.0027*1.75,0.0004*1.75,6.12e3*1.75,2.3e3,2.3e3]);%1.85e-4;%3e-4;% 
end

%--------------------------%
trace_transit =0; % Variable for dranwing the transition poses
Const.q = q;
nb_modes_theta = Const.dim_base_q_theta/(Config.nb_tubes-1);


q_e = Const.q(1:Const.dim_base_q_e);
q_t = Const.q(Const.dim_base_q_e+1:end);

Const.q_e = q_e;
Const.q_t = q_t;
N_noeuds = Config.N_noeuds;
Const.it_troncon = 1;

[DX,X_grille]=cheb(N_noeuds-1,Config.Li(1));  % sur une grille [0,L]
X_k{1} = X_grille;

% Calculation of X_k on the X1 grid

for it_tubes = 2 : Config.nb_tubes 
    Config.it_tubes = it_tubes;
    [T X_k{it_tubes}] = ode45(@(X_1,X_k) face_2_face(X_1,X_k,Const,Config),flip(X_grille),0); % Solve ODE
    X_k{it_tubes} = flip(X_k{it_tubes});
end

% Calculation of Xi_1 on the X1 grid
for it_x = 1:N_noeuds
    Phi_e(:,:,it_x) = Base_Phi(X_grille(it_x),time,Const,Config);
    
    Xi_K = Phi_e(:,:,it_x)'*Const.q(1:Const.dim_base_q_e);
    Xi_a = [Xi_K;0;0;Config.D{1}(X_grille(it_x))*Xi_K(1)];
    Xi_0(:,it_x) = [Config.K_0{1}(X_grille(it_x));0;0;1;Config.D_prime{1}(X_grille(it_x));Config.D{1}(X_grille(it_x))*Config.K_0{1}(X_grille(it_x))];
    Xi(:,it_x) = Xi_a +  Xi_0(:,it_x);
    for it_tube = 2:Config.nb_tubes
        [Phi_theta(:,:,it_x),Phi_theta_prime(:,:,it_x)] = Base_Phi_prime(X_k{it_tubes}(it_x),time,Const,Config);
        [Phi_theta(:,:,it_x),Phi_theta_prime(:,:,it_x)] = Base_Phi_prime(Config.L-X_grille(it_x),time,Const,Config);
        Delta_theta{it_tube}(:,it_x) = Phi_theta(:,:,it_x)'*q_t(1+(it_tube-2)*nb_modes_theta:(it_tube-1)*nb_modes_theta);
        Delta_theta_prime{it_tube}(:,it_x) = Phi_theta_prime(:,:,it_x)'*q_t(1+(it_tube-2)*nb_modes_theta:(it_tube-1)*nb_modes_theta);
        h_k{it_tube}(it_x) =  face_2_face(X_k{1}(it_x),X_k{it_tubes}(it_x),Const,Config);
        h_k_0{it_tube}(it_x) =  face_2_face_0(X_k{1}(it_x),X_k{it_tubes}(it_x),Const,Config);
    end
end

Xi_k{1} = Xi;
Xi_k_0{1} = Xi_0;
Const.Xi = Xi;

Delta_theta{1} = 0*Delta_theta{2};
Delta_theta_prime{1} = 0*Delta_theta_prime{2};


% Quaternion of the beam 1

CI = q_0;
CL_ind=[1,N_noeuds+1,2*N_noeuds+1,3*N_noeuds+1];

f_A = @(it) quaternion_dot2(Xi(:,it));
f_B = @(it) zeros(max(size(CI)),1);

QX_k{1} = integral_spectral(f_A,f_B,CI,DX,N_noeuds,CL_ind);

q_0 = QX_k{1}(:,end);
q_0_1 = q_0;

%  position of the beam 1

CI = r_0;
CL_ind=[1,N_noeuds+1,2*N_noeuds+1,3*N_noeuds+1];

f_A = @(it) zeros(max(size(CI)),max(size(CI)));
f_B = @(it) r_dot(QX_k{1}(:,it),Xi(4:6,it));

rX{1} = integral_spectral(f_A,f_B,CI,DX,N_noeuds,CL_ind);

r_0 = rX{1}(:,end);

if trace_transit==1
    figure(1)
    plot3(rX{1}(1,:),rX{1}(2,:),rX{1}(3,:),'LineWidth',2)
    hold on
    grid on
    axis equal
    xlabel('X (m)')
    ylabel('Y (m)')
    zlabel('Z (m)')
end


for it_tubes = 2 : Config.nb_tubes 
    
    % Calculation of Xi_k on the X1 grid
    
    for it_x = 1:N_noeuds
        theta_k = Config.theta_0{it_tubes}(X_k{it_tubes}(it_x))+Delta_theta{it_tube}(:,it_x);
        Rx = [1 0 0;0 cos(theta_k) -sin(theta_k);0 sin(theta_k) cos(theta_k)];
        Rz = [-1, 0, 0;0, -1, 0;0, 0, 1];
        r1_k = Rz*[((Config.D{1}(X_grille(it_x))*eye(3)+Config.D{it_tubes}(X_k{it_tubes}(it_x))*Rx)*[0;1;0])];
        g1_k = [Rz*Rx,r1_k; 0 0 0 1];
        gk_1 = inv(g1_k);
        Rk_1 = gk_1(1:3,1:3);
        rk_1 = gk_1(1:3,4);
        
        mu_k = -[1 0 0 0 0 Config.D{it_tubes}(X_k{it_tubes}(it_x))]';
        theta_prime_0_k = Config.theta_0_prime{it_tubes}(X_k{it_tubes}(it_x));
        
        nu_k = [theta_prime_0_k,0,0,0,cos(theta_k)*Config.D_prime{1}(X_grille(it_x))-Config.D_prime{it_tubes}(X_k{it_tubes}(it_x)),sin(theta_k)*Config.D_prime{1}(X_grille(it_x))+Config.D{it_tubes}(X_k{it_tubes}(it_x))*theta_prime_0_k]';
        
        Xi_k_1 = mu_k * Delta_theta_prime{it_tube}(:,it_x) + nu_k;
        
        Xi_k{it_tubes}(:,it_x) = (Ad_g_(Rk_1,rk_1)*Xi(:,it_x) + Xi_k_1); 
    end
    
    % Quaternion of the beam k
    theta = Config.theta_0{it_tubes}(X_k(end))+Delta_theta{it_tube}(:,end);
    Rx = [1 0 0;0 cos(theta_k) -sin(theta_k);0 sin(theta_k) cos(theta_k)];
    Rz = [-1, 0, 0;0, -1, 0;0, 0, 1];
    r1_k = Rz*[((Config.D{1}(X_grille(end))*eye(3)+Config.D{it_tubes}(X_k{it_tubes}(end))*Rx)*[0;1;0])];
    g1_k = [Rz*Rx,r1_k; 0 0 0 1];
    g_0_1 = [quat2rot(q_0_1),r_0;0 0 0 1];
    g_0_k = g_0_1*g1_k;
    
    CI = rot2quat(g_0_k(1:3,1:3));
    CL_ind=[N_noeuds,2*N_noeuds,3*N_noeuds,4*N_noeuds];

    f_A = @(it) quaternion_dot2(Xi_k{it_tubes}(:,it));
    f_B = @(it) zeros(max(size(CI)),1);

    QX_k{it_tubes} = integral_spectral(f_A,f_B,CI,DX,N_noeuds,CL_ind);

    % Position of the beam 2
       
    CI = g_0_k(1:3,4);
    CL_ind=[N_noeuds,2*N_noeuds,3*N_noeuds];

    f_A = @(it) zeros(max(size(CI)),max(size(CI)));
    f_B = @(it) r_dot(QX_k{it_tubes}(:,it),Xi_k{it_tubes}(4:6,it));
    
    rX{it_tubes} = integral_spectral(f_A,f_B,CI,DX,N_noeuds,CL_ind);

    am(it_tubes)=Config.L-X_k{it_tubes}(1); % a
    thetam(it_tubes) = Delta_theta{it_tube}(:,1);

    if trace_transit==1
        figure(1)
        plot3(rX{it_tubes}(1,:),rX{it_tubes}(2,:),rX{it_tubes}(3,:),'LineWidth',2)
    end


end


for it_tubes = Config.nb_tubes :-1:  2
    % Calculation of Lambda of the beam 2 X1: 0-> l
    Const.it_troncon = 1;
    for it_x = 1:N_noeuds
        [~,F_1_e,F_bar(:,it_x)] = Forces_exterieur(X_grille(it_x),QX_k{it_tubes}(:,it_x),rX{it_tubes}(:,it_x),0,Const,Config);
    end

    F_1 = [Const.C{it_tubes}; 0; 0; Const.T{it_tubes}; 0; 0];
    CI = F_1;

    CL_ind=[1,N_noeuds+1,2*N_noeuds+1,3*N_noeuds+1,4*N_noeuds+1,5*N_noeuds+1];

    f_A = @(it) ad_(Xi_k{it_tubes}(:,it))';
    f_B = @(it) -h_k{it_tubes}(it)*F_bar(:,it);
    
    Lambda_X_k{it_tubes} = integral_spectral(f_A,f_B,CI,DX,N_noeuds,CL_ind);
    F_k_p{it_tubes} = Lambda_X_k{it_tubes}(:,end);

end


    % Calculation of Lambda of the beam 1 X1: 1-> 0
Const.it_troncon = 1;
for it_x = 1:N_noeuds
    [~,F_1_e,F_bar(:,it_x)] = Forces_exterieur(X_grille(it_x),QX_k{1}(:,it_x),rX{1}(:,it_x),0,Const,Config);
end
F_1_p = F_1_e;

    for it_tubes = Config.nb_tubes :-1:  2
        theta = Config.theta_0{it_tubes}(X_k(end))+Delta_theta{it_tube}(:,end);
        Rx = [1 0 0;0 cos(theta_k) -sin(theta_k);0 sin(theta_k) cos(theta_k)];
        Rz = [-1, 0, 0;0, -1, 0;0, 0, 1];
        r1_k = Rz*[((Config.D{1}(X_grille(end))*eye(3)+Config.D{it_tubes}(X_k{it_tubes}(end))*Rx)*[0;1;0])];
        g1_k = [Rz*Rx,r1_k; 0 0 0 1];
        R1_k = Rz*Rx;
        F_1_p = inv(Ad_g_(R1_k,r1_k))'*F_k_p{it_tubes}(:,end)+F_1_p;
    end

CI = F_1_p;

CL_ind=[N_noeuds,2*N_noeuds,3*N_noeuds,4*N_noeuds,5*N_noeuds,6*N_noeuds];

f_A = @(it) ad_(Xi(:,it))';
f_B = @(it) -F_bar(:,it);

Lambda_X_k{1} = integral_spectral(f_A,f_B,CI,DX,N_noeuds,CL_ind);
F_k_p{1} = Lambda_X_k{1}(:,1);
 
%% ------------------------------------%

% Calculation of Qa and residu  X1: 0-> l

[Xi_k_0] = function_Xi_0(time,Const,Config);
    
for it_tubes = 1 : Config.nb_tubes     
    for it_x = 1:N_noeuds
        
        hk =  face_2_face(X_k{1}(it_x),X_k{it_tubes}(it_x),Const,Config);
        hk_0 =  face_2_face_0(X_k{1}(it_x),X_k{it_tubes}(it_x),Const,Config);

        if it_tubes==1
            hk = 1;
            hk_0 = 1;
        end
        Lambda_X_k{it_tubes}(:,it_x) =  Lambda_X_k{it_tubes}(:,it_x) - Hook{it_tubes}(X_k{it_tubes}(it_x))*((1/hk)*Xi_k{it_tubes}(:,it_x) - (1/hk_0)*Xi_k_0{it_tubes}(:,it_x)); 
        
        theta_k = Config.theta_0{it_tubes}(X_k{it_tubes}(it_x))+Delta_theta{it_tube}(:,it_x);
        Rx = [1 0 0;0 cos(theta_k) -sin(theta_k);0 sin(theta_k) cos(theta_k)];
        Rz = [-1, 0, 0;0, -1, 0;0, 0, 1];
        r1_k = Rz*[((Config.D{1}(X_grille(it_x))*eye(3)+Config.D{it_tubes}(X_k{it_tubes}(it_x))*Rx)*[0;1;0])];
        g1_k = [Rz*Rx,r1_k; 0 0 0 1];
        gk_1 = inv(g1_k);
        Rk_1 = gk_1(1:3,1:3);
        rk_1 = gk_1(1:3,4);
        
        delta_h_k_xi_1 = -[0, sin(theta)*Config.D{it_tubes}(X_k{it_tubes}(it_x)),(Config.D{1}(X_grille(it_x))+cos(theta)*Config.D{it_tubes}(X_k{it_tubes}(it_x))),0,0,0]';
         
        coef_dual_epsilon{it_tubes}(:,:,it_x) = Ad_g_(Rk_1,rk_1)' - delta_h_k_xi_1*(1/hk)*Xi_k{it_tubes}(:,it_x)';   

        [Phi,Phi_prime] = Base_Phi_prime(Config.L-X_grille(it_x),time,Const,Config);;
        nb_modes_theta = Const.dim_base_q_theta/(Config.nb_tubes-1);
        mu_k = -[1 0 0 0 0 Config.D{it_tubes}(X_k{it_tubes}(it_x))]';
        
        coef_dual_theta_prime = Phi_prime*mu_k'*Lambda_X_k{it_tubes}(:,it_x);
        theta_0_k = Config.theta_0{it_tubes}(X_k{it_tubes}(it_x));
        delta_h_delta_theta_k = Config.D{it_tubes}(X_grille(it_x)) * (Xi(3,it_x)*sin(theta_0_k) - Xi(2,it_x)*cos(theta_0_k));
        delta_nu_k_delta_theta_k = Config.D_prime{1}(X_grille(it_x))*[0;0;0;0;-sin(theta_0_k);cos(theta_0_k)]; 
        coef_dual_theta = Phi *((delta_nu_k_delta_theta_k' + mu_k'*ad_(Ad_g_(Rk_1,rk_1)'*Xi(:,it_x)) - delta_h_delta_theta_k*Xi_k{it_tubes}(:,it_x)')*Lambda_X_k{it_tubes}(:,it_x));
        Lambda_X_theta{it_tubes}(:,it_x)=(coef_dual_theta_prime+coef_dual_theta);

    end
    
end

Total_Lambda_X_k = 0*Lambda_X_k{1};
Total_Lambda_X_theta = [];
for it_tubes = 2 : Config.nb_tubes 
    for it_x = 1:N_noeuds
        Total_Lambda_X_k(:,it_x) = Total_Lambda_X_k(:,it_x)+coef_dual_epsilon{it_tubes}(:,:,it_x)*Lambda_X_k{it_tubes}(:,it_x);
    end
    Total_Lambda_X_theta = [Total_Lambda_X_theta;Lambda_X_theta{it_tubes}];
end
    CI =  zeros(Const.dim_base_q_e,1);
    CL_ind=[];
    for i=1:Const.dim_base_q_e
        CL_ind=[CL_ind,(i-1)*N_noeuds+1];
    end

    f_A = @(it) zeros(Const.dim_base_q_e,Const.dim_base_q_e);
    f_B = @(it) [Phi_e(:,:,it),zeros(Const.dim_base_q_e,2),Phi_e(:,1,it)*Config.D{1}(X_grille(it))]*(-Lambda_X_k{1}(:,it) + Total_Lambda_X_k(:,it));

    Qa_k_X = integral_spectral(f_A,f_B,CI,DX,N_noeuds,CL_ind);

    Qa_T_e = Qa_k_X(:,end);
    
    
% Calculation of Qtheta and residu  X1: 0-> l

CI =  zeros(Const.dim_base_q_theta,1);
CL_ind=[];
for i=1:Const.dim_base_q_theta
    CL_ind=[CL_ind,(i-1)*N_noeuds+1];
end

f_A = @(it) zeros(Const.dim_base_q_theta,Const.dim_base_q_theta);
f_B = @(it)  Total_Lambda_X_theta(:,it);

Qa_t_X = integral_spectral(f_A,f_B,CI,DX,N_noeuds,CL_ind);

Qa_T_t = Qa_t_X(:,end);

Q_a_T = [Qa_T_e;-Qa_T_t];
