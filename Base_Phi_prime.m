function [Phi,Phi_prime] = Base_Phi_prime(X,t,Const,Config)
    
    
    X = X/Config.Li(Const.it_troncon);

    Base = [2*X;-6*X+6*X^2;12*X-30*X^2+30*X^3]/2;
    Base_prime = [2;-6+12*X;12-60*X+90*X^3]/(2*Config.Li(Const.it_troncon));
    Phi = [];
    Phi_prime = [];
     
    Phi = Base(1:Const.dim_base_q_theta/(Config.nb_tubes-1));
    Phi_prime = Base_prime(1:Const.dim_base_q_theta/(Config.nb_tubes-1));
