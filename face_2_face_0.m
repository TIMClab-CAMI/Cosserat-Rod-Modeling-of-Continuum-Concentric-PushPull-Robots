function  X_k_prime = face_2_face_0(X_1,X_k,Const,Config)

it_tube = Config.it_tubes;

q_e = 0*Const.q(1:Const.dim_base_q_e);
q_t = 0*Const.q(Const.dim_base_q_e+1:end);

D1 = Config.D{1};
D_k = Config.D{it_tube};

nb_modes_theta = Const.dim_base_q_theta/(Config.nb_tubes-1);

[Phi,~] = Base_Phi_prime(Config.L-X_1,0,Const,Config);
Delta_theta_k = Phi'*q_t(1+(it_tube-2)*nb_modes_theta:(it_tube-1)*nb_modes_theta);

Theta_k = Config.theta_0{it_tube}(X_k)+Delta_theta_k;

Xi = Const.B*(Base_Phi(X_1,0,Const,Config)'*q_e) +  [Config.K_0{1}(X_1);0;0;1;Config.D_prime{1}(X_1);Config.D{1}(X_1)*Config.K_0{1}(X_1)];

X_k_prime = -(1 + D1(X_1)* Xi(3) + D_k(X_k)*(cos(Theta_k)* Xi(3) + sin(Theta_k)* Xi(2)));
