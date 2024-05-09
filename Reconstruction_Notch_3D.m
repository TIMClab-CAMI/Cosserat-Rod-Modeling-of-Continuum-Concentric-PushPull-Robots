function   [QX_T,rX_T] = Reconstruction_Notch_3D(q_0,r_0,Const,Config)
time =0;

D = Config.D;
D_prime = Config.D_prime;
N_nodes =200; 
X_grid_T = [];
Xi_T = [];
time = 0;

trace_transit =0;

%--------------------------%
nb_modes_theta = Const.dim_base_q_theta/(Config.nb_tubes-1);

q_e = Const.q(1:Const.dim_base_q_e);
q_t = Const.q(Const.dim_base_q_e+1:end);

Const.q_e = q_e;
Const.q_t = q_t;
Const.it_troncon = 1;

[DX,X_grid]=cheb(N_nodes-1,Config.Li(1));  % sur une grille [0,L]
X_k{1} = X_grid;

% Calculation of X_k in X1

for it_tubes = 2 : Config.nb_tubes 
    Config.it_tubes = it_tubes;
    [T X_k{it_tubes}] = ode45(@(X_1,X_k) face_2_face(X_1,X_k,Const,Config),flip(X_grid),0); % Solve ODE
    X_k{it_tubes} = flip(X_k{it_tubes});
end

% Calculation of Xi_1 in X_1 of beam 1 
for it_x = 1:N_nodes
    Phi_e(:,:,it_x) = Base_Phi(X_grid(it_x),time,Const,Config);
    
    Xi_K = Phi_e(:,:,it_x)'*Const.q(1:Const.dim_base_q_e);
    Xi_a = [Xi_K;0;0;Config.D{1}(X_grid(it_x))*Xi_K(1)];
    Xi_0(:,it_x) = [Config.K_0{1}(X_grid(it_x));0;0;1;Config.D_prime{1}(X_grid(it_x));Config.D{1}(X_grid(it_x))*Config.K_0{1}(X_grid(it_x))];
    Xi(:,it_x) = Xi_a +  Xi_0(:,it_x);
    for it_tube = 2:Config.nb_tubes
        [Phi_theta(:,:,it_x),Phi_theta_prime(:,:,it_x)] = Base_Phi_prime(X_k{it_tubes}(it_x),time,Const,Config);
        [Phi_theta(:,:,it_x),Phi_theta_prime(:,:,it_x)] = Base_Phi_prime(Config.L-X_grid(it_x),time,Const,Config);
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


% Quaternion of beam 1

CI = q_0;
CL_ind=[1,N_nodes+1,2*N_nodes+1,3*N_nodes+1];

f_A = @(it) quaternion_dot2(Xi(:,it));
f_B = @(it) zeros(max(size(CI)),1);

QX_k{1} = integral_spectral(f_A,f_B,CI,DX,N_nodes,CL_ind);

q_0 = QX_k{1}(:,end);
q_0_1 = q_0;

% Position of beam 1

CI = r_0;
CL_ind=[1,N_nodes+1,2*N_nodes+1,3*N_nodes+1];

f_A = @(it) zeros(max(size(CI)),max(size(CI)));
f_B = @(it) r_dot(QX_k{1}(:,it),Xi(4:6,it));

rX{1} = integral_spectral(f_A,f_B,CI,DX,N_nodes,CL_ind);

r_0 = rX{1}(:,end);

if trace_transit==1
    h=figure;
    plot3(-rX{1}(2,:),-rX{1}(3,:),rX{1}(1,:),'LineWidth',2)
    hold on
    grid on
    axis equal
    xlabel('X (m)')
    ylabel('Y (m)')
    zlabel('Z (m)')

end


for it_tubes = 2 : Config.nb_tubes 
    % Calculation of Xi_k in X1
    
    for it_x = 1:N_nodes
        theta_k = Config.theta_0{it_tubes}(X_k{it_tubes}(it_x))+Delta_theta{it_tube}(:,it_x);
        Rx = [1 0 0;0 cos(theta_k) -sin(theta_k);0 sin(theta_k) cos(theta_k)];
        Rz = [-1, 0, 0;0, -1, 0;0, 0, 1];
        r1_k = Rz*[((Config.D{1}(X_grid(it_x))*eye(3)+Config.D{it_tubes}(X_k{it_tubes}(it_x))*Rx)*[0;1;0])];
        g1_k = [Rz*Rx,r1_k; 0 0 0 1];
        gk_1 = inv(g1_k);
        Rk_1 = gk_1(1:3,1:3);
        rk_1 = gk_1(1:3,4);
        
        mu_k = -[1 0 0 0 0 Config.D{it_tubes}(X_k{it_tubes}(it_x))]';
        theta_prime_0_k = Config.theta_0_prime{it_tubes}(X_k{it_tubes}(it_x));
        
        nu_k = [theta_prime_0_k,0,0,0,cos(theta_k)*Config.D_prime{1}(X_grid(it_x))-Config.D_prime{it_tubes}(X_k{it_tubes}(it_x)),sin(theta_k)*Config.D_prime{1}(X_grid(it_x))+Config.D{it_tubes}(X_k{it_tubes}(it_x))*theta_prime_0_k]';
        
        Xi_k_1 = mu_k * Delta_theta_prime{it_tube}(:,it_x) + nu_k;
        
        Xi_k{it_tubes}(:,it_x) = (1/h_k{it_tube}(it_x))*(Ad_g_(Rk_1,rk_1)*Xi(:,it_x) + Xi_k_1); 
    end
    
    % Quaternion of beam k
    theta = Config.theta_0{it_tubes}(X_k(end))+Delta_theta{it_tube}(:,end);
    Rx = [1 0 0;0 cos(theta_k) -sin(theta_k);0 sin(theta_k) cos(theta_k)];
    Rz = [-1, 0, 0;0, -1, 0;0, 0, 1];
    r1_k = Rz*[((Config.D{1}(X_grid(end))*eye(3)+Config.D{it_tubes}(X_k{it_tubes}(end))*Rx)*[0;1;0])];
    g1_k = [Rz*Rx,r1_k; 0 0 0 1];
    g_0_1 = [quat2rot(q_0_1),r_0;0 0 0 1];
    g_0_k = g_0_1*g1_k;
    
    CI = rot2quat(g_0_k(1:3,1:3));
    CL_ind=[N_nodes,2*N_nodes,3*N_nodes,4*N_nodes];

    f_A = @(it) quaternion_dot2(h_k{it_tube}(it)*Xi_k{it_tubes}(:,it));
    f_B = @(it) zeros(max(size(CI)),1);

    QX_k{it_tubes} = integral_spectral(f_A,f_B,CI,DX,N_nodes,CL_ind);

    % Position of beam k
       
    CI = g_0_k(1:3,4);
    CL_ind=[N_nodes,2*N_nodes,3*N_nodes];

    f_A = @(it) zeros(max(size(CI)),max(size(CI)));
    f_B = @(it) r_dot(QX_k{it_tubes}(:,it),h_k{it_tube}(it)*Xi_k{it_tubes}(4:6,it));
    
    rX{it_tubes} = integral_spectral(f_A,f_B,CI,DX,N_nodes,CL_ind);

    am(it_tubes)=Config.L-X_k{it_tubes}(1); % a
    thetam(it_tubes) = Delta_theta{it_tube}(:,1);

    if trace_transit==1
        plot3(-rX{it_tubes}(2,:),-rX{it_tubes}(3,:),rX{it_tubes}(1,:),'LineWidth',2)
    end
    alpha_trace = [-pi:pi/20:pi];
    it_tubes = 2;
end

if trace_transit==1
    hold off
    axis equal
    xlabel('X (m)')
    ylabel('Y (m)')
    zlabel('Z (m)')
    legend('tube 1','tube 2')
end
QX_T = QX_k;
rX_T = rX;
drawnow
%-----------------------%

L =Config.L;

X = X_grid;
pas_dent = 0.005;
for it=1:max(size(X))
    theta = [-pi:pi/20:pi];

    if (mod(round(X_k{1}(it)/pas_dent),2)==0)
        theta_max = pi;
    else 
        theta_max = pi/4;
    end
    ind = find((abs(theta)<=theta_max)==0);
    theta(ind) = nan;
    X_1_0 = 0*X(it)*ones(size(theta));
    Y_1_0 = 3e-3*cos(theta)-Config.D{1}(rX{1}(:,it));
    Z_1_0 = 3e-3*sin(theta);
    
    g_0_1 = [quat2rot(QX_k{1}(:,it)),rX{1}(:,it);0 0 0 1];
    for it_theta = 1 :max(size(theta))
        temp = g_0_1*[X_1_0(it_theta );Y_1_0(it_theta );Z_1_0(it_theta );1];
        X_1(it,it_theta ) = temp(1);
        Y_1(it,it_theta ) = temp(2);
        Z_1(it,it_theta ) = temp(3);
    end
    theta = [-pi:pi/20:pi];

    if (mod(round(X_k{2}(it)/pas_dent),2)==0)
        theta_max = pi;
    else 
        theta_max = pi/4;
    end
    ind = find((abs(theta)<=theta_max)==0);
    theta(ind) = nan;
    X_2_0 = 0*X(it)*ones(size(theta));
    Y_2_0 = 1.5e-3*cos(theta)-Config.D{2}(rX{2}(:,it));
    Z_2_0 = 1.5e-3*sin(theta);
    
    g_0_2 = [quat2rot(QX_k{2}(:,it)),rX{2}(:,it);0 0 0 1];
    for it_theta = 1 :max(size(theta))
        temp = g_0_2*[X_2_0(it_theta );Y_2_0(it_theta );Z_2_0(it_theta );1];
        X_2(it,it_theta ) = temp(1);
        Y_2(it,it_theta ) = temp(2);
        Z_2(it,it_theta ) = temp(3);
    end
end
figure(200)
surf(-Y_1,-Z_1,X_1,'EdgeColor','none','FaceAlpha',0.5,'Facecolor',[0.41,0.81,0.70])
hold on
surf(-Y_2,-Z_2,X_2,'EdgeColor','none','FaceAlpha',1,'Facecolor',[0.96,0.40,0.33])
grid on 

surf(0.06*ones(size(Y_1)),-Z_1,X_1,'EdgeColor','none','FaceAlpha',1,'Facecolor',0.8*[1,1,1])
surf(-Y_1,0.06*ones(size(Z_1)),X_1,'EdgeColor','none','FaceAlpha',1,'Facecolor',0.8*[1,1,1])
surf(-Y_1,-Z_1,0.*ones(size(Z_1)),'EdgeColor','none','FaceAlpha',1,'Facecolor',0.8*[1,1,1])

surf(0.06*ones(size(Y_2)),-Z_2,X_2,'EdgeColor','none','FaceAlpha',1,'Facecolor',0.8*[1,1,1])
surf(-Y_2,0.06*ones(size(Z_2)),X_2,'EdgeColor','none','FaceAlpha',1,'Facecolor',0.8*[1,1,1])
surf(-Y_2,-Z_2,0.*ones(size(Z_2)),'EdgeColor','none','FaceAlpha',1,'Facecolor',0.8*[1,1,1])
hold off
axis equal
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
xlim([-0.06 0.06])
ylim([-0.06 0.06])
zlim([0 0.160])

legend('tube 1','tube 2','Orientation','horizontal','Location','north')
drawnow

 