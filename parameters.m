
Const.rho = 1000;
Const.Rc  = 3.5e-3;

Aire = pi*Const.Rc^2;
Const.Aire = Aire;

% --------------------- %

V_a = Config.V_a;

M_selec = eye(6,6);
M_selec_K = zeros(3,3);

[row,col] = find(V_a==1);

Const.B = M_selec(:,col);
Const.B_bar = M_selec;
Const.B_bar(:,col)=[];

[row,col] = size(col);

Const.dim_B = col;

% Imposed strain
Const.Xi_c = Const.B_bar'*[0;0;0;1;0;0];
Const.Xi_0 = Const.B'*[0;0;0;1;0;0];


