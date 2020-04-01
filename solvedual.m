function [nu_star, H_star, Z_star] = solvedual(M0, K)
%% description
% solves max nu 
%        st  M0 - K'HK \geq nu I
%        H doubly hyperdominant
% as expressed in Section IVb of 
%   1) Lee, Seiler Finite Step Performance Guarantees for First Order Optimization
%   Algorithms: Revisited from a Control Theoretic Perspective, arxiv 2020
%% Yalmip
N = size(K,1)/2;
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

constraints = [constraints M0-K'*M*K >= nu*eye(n)];

%maximize nu by minimizing -nu
ops = sdpsettings('solver','mosek','verbose',0,'debug',0);
optimize(constraints, -nu, ops);

nu_star = value(nu);
H_star = value(H);
M_star = value(M);
Z_star = M0-K'*M_star*K-nu_star*eye(n);