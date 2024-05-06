function [T,q] = Time_Integration_Newton_Static_Notch_Spectral(t,Const,Config)

% Newton residual criterion
r_min = 1e-10;

% Output Initialisation
q(:,1)         = Const.q;
T(1) = 0;

% Position and orietation of the Base of CCPPR 
q_0   = Const.q_0;
r_0   = Const.r_0;

% Définition of the desired control values
a_2 = [-10:1:10];
a_3 = a_2;
thetam_2 = 0;
thetam_3 = thetam_2;

% Number of simulation steps
N_simu = max(size(a_2));

for it_simu = 1:N_simu
    for it_theta = 1:max(size(thetam_2))
    time = 0;
    for it_tubes = 2 : Config.nb_tubes  
        Const.T{it_tubes} = 0;
        Const.C{it_tubes} = 0;
    end
    
    % Assignment of the current state for the Newton loop at t+1 (np1 => n+1)
    q_n = 0*q(:,it_simu);
    % State prediction
    q_np1_k = q_n;
    
    % Assignation of the desired values
    ad(2) =  a_2(it_simu)*1e-3;
    ad(3) =  a_3(it_simu)*1e-3;
    thetam_d(2) =  thetam_2(it_theta);
    thetam_d(3) =  thetam_3(it_theta);
    
    % Initialisation of the measured value
    am = 10*ones(1,Config.nb_tubes);
    thetam = 10*ones(1,Config.nb_tubes);
    
    for it_tubes = 2 : Config.nb_tubes  
        e(it_tubes) = (am(it_tubes)-ad(it_tubes));
        e_theta(it_tubes) = (thetam(it_tubes)-thetam_d(it_tubes));
    end
    
    while max(abs([e,e_theta*1e-2]))>1e-4

    Const.q  = q_np1_k;

    % ISM output and Jacobian at t+1 
    [Q_a_T] = ISM_spectral_notch(q_0,r_0,Const.q,0,0,time,Const,Config);
    r = Q_a_T;

    % root solution     

    [q_np1_k, Feval,EXITFLAG] =  fsolve(@(x) ISM_spectral_notch(q_0,r_0,x,0,0,time,Const,Config), q_np1_k,optimoptions('fsolve','Display','iter-detailed','Jacobian','off','TolFun',1e-12,'StepTolerance',1e-8,'Algorithm','levenberg-marquardt'));%,'InitDamping',0))
    Const.q  = q_np1_k;
 
   [~,am,thetam] = ISM_spectral_notch_a(q_0,r_0,Const.q,0,0,time,Const,Config);

   % PID Controleur 
    for it_tubes = 2 : Config.nb_tubes  
        Const.T{it_tubes} = Const.T{it_tubes} - 2e2*(am(it_tubes)-ad(it_tubes));
        Const.C{it_tubes} = Const.C{it_tubes} - 1e-2*(thetam(it_tubes)-thetam_d(it_tubes));
    end
    
    q_np1_kp1 = q_np1_k;
    for it_tubes = 2 : Config.nb_tubes  
        e(it_tubes) = (am(it_tubes)-ad(it_tubes));
        e_theta(it_tubes) = (thetam(it_tubes)-thetam_d(it_tubes));
    end

    end
    % Drawing
    [QX_T,rX_T] = Reconstruction_Notch_3D(q_0,r_0,Const,Config);
    title(['a = ' num2str(ad(2)*1000) 'mm, \theta = ' num2str(thetam_d(2) *180/pi) '°'])

    % Assignment of the converged state when r<r_min
    q(:,it_simu+1)         = q_np1_kp1;
    T(it_simu+1) = time;

    end
end
