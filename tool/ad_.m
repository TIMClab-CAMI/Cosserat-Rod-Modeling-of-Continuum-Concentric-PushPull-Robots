function M=ad_(G)

W=G(1:3);
U=G(4:6);

ad11=hat_(W);
ad12=zeros(3);
ad21=hat_(U);

M=[[ad11,ad12];[ad21,ad11]];
