function [nu_star, H_star, Z_star] = solvedual(M0, R2)
%% description
% solves max nu 
%        st  M0 - R2' [0 H; H' 0] R2 \geq nu I
%        H doubly hyperdominant
%% Yalmip
N = size(R2,1)/2;
n = size(M0,1);
nu = sdpvar(1);
H = sdpvar(N,N,'full');
M = [zeros(N) H; H' zeros(N)];

constraints = [];

%constrain H to be doubly hyperdominant
%off diagonal entries nonpositive
for i=1:N
    for j=1:N
        if i~=j
            constraints = [constraints H(i,j)<=0];
        end
    end
end

%row and column sumns nonnegative
for i=1:N
    constraints = [constraints sum(H(i,:))>=0];
    constraints = [constraints sum(H(:,i))>=0];
end

constraints = [constraints M0-R2'*M*R2 >= nu*eye(n)];

%maximize nu by minimizing -nu
ops = sdpsettings('solver','mosek','verbose',0,'debug',0);
optimize(constraints, -nu, ops);

nu_star = value(nu);
H_star = value(H);
M_star = value(M);
Z_star = M0-R2'*M_star*R2-nu_star*eye(n);