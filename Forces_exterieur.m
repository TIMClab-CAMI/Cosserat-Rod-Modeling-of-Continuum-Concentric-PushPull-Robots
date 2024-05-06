function [F_0,F_1,F_bar] = Forces_exterieur(X,q,r,eta,Const,Config)
it_troncon = Const.it_troncon;

R = quat2rot(q);
g = [R,r;0 0 0 1];

F_bar = Const.Fbar_materielle + [R',zeros(3,3);zeros(3,3),R']*Const.Fbar_spaciale;
F_0 = zeros(6,1);
F_1 = zeros(6,1);

if X <= 0
    F_0   = Const.Fm_materielle + [R',zeros(3,3);zeros(3,3),R']*Const.Fm_spaciale;
    F_bar = zeros(6,1);
end
if X >= Config.Li(it_troncon)
    r_pully = [0.150;-0.024;0];
    beta = atan2(r(2)-r_pully(2),r(1)-r_pully(1));
    Fp_spaciale = 0*[0;0;0;20e-3*9.81*[cos(beta);sin(beta);0]];
    Fp_spaciale = 0*[0;0;0;0;-50e-3*9.81;0];
    Fp_spaciale = 0*[0;0;0;0;0;50e-3*9.81];
    F_1   = Const.Fp_materielle + [R',zeros(3,3);zeros(3,3),R']*Fp_spaciale;
    F_bar = zeros(6,1);
end


