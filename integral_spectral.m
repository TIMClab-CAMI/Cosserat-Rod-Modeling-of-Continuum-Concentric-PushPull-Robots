function X = integral_spectral(f_A,f_B,CI,DX,N_noeuds,CL_ind)

n = max(size(CI));

A = zeros(n*N_noeuds,n*N_noeuds);
b = zeros(n*N_noeuds,1);

D = [];
for i=1:n
    D=blkdiag(D,DX);
end

P=eye(size(A));
for i=1:n
    P_sauv=P(:,i);
    P(:,i)=P(:,CL_ind(i));
    P(:,CL_ind(i))=P_sauv;
end
    
for it_x = 1:N_noeuds
    
    A_xi  = f_A(it_x);
    b_xi = f_B(it_x);

    for it_l = 1:n
        for it_c = 1:n
            A(it_x + (it_l-1)*N_noeuds,it_x + (it_c-1)*N_noeuds) = A_xi(it_l,it_c);
        end
        b(it_x + (it_l-1)*N_noeuds,1) = b_xi(it_l,1);
    end
end

A=P'*A*P;
b=P'*b;
D=P'*D*P;

CL  = (D(:,1:n)-A(:,1:n))*CI;
res = (D(n+1:end,n+1:end)-A(n+1:end,n+1:end))\(b(n+1:end,1)-CL(n+1:end,1));
X_x=P*[CI;res];
X_x_M = reshape(X_x,[N_noeuds,n]);
X = X_x_M';