% Solves bound from Theorem 2 in 
%   1) Lee, Seiler. Finite Step Performance of First-order Methods Using Interpolation
%   Conditions Without Function Evaluations, arxiv 2020
%  For heavy ball optimized for quadratics.

%% Define Problem Class
%strongly convex with parameter m, L-Lipschitz with parameter L
% Assumes that L,m and T are defined in workspace
k = L/m;

%% Define algorithm
%heavy ball optimized for quadratics

alpha = 4/(sqrt(L)+sqrt(m))^2;
beta = ((sqrt(k)-1)/(sqrt(k)+1))^2;

A = [1+beta,-beta;1,0];
B = [-alpha;0];
C = [1,0];
D = [0];
%% Define IQC filter
%cyclic pointwise IQC
S = [L -1; -m 1];
N = T;

%N is number of (u,y) pairs
%[tilde{y}; tilde{u}] = S [y; u]
%want z_k = [tilde{y}_k; tilde{y}_{k-1}; ...; tilde{y}_{k-N};
%               tilde{u}_k; tilde{u}_{k-1}; ...; tilde{u}_{k-N}]
if N==0
    Apsi = [];
    Bpsi = zeros(0,2);
else    
    Apsi = [zeros(1, 2*N); eye(N-1) zeros(N-1, N+1); zeros(1,2*N); zeros(N-1,N) eye(N-1), zeros(N-1, 1)];
    Bpsi = [S(1,:); zeros(N-1,2); S(2,:); zeros(N-1,2)];
end

Cpsi = [zeros(1, 2*N); eye(N) zeros(N, N); zeros(1, 2*N); zeros(N,N) eye(N)];
Dpsi = [S(1,:); zeros(N, 2); S(2,:); zeros(N,2)];
%% Define Augmented system
Bpsi_y = Bpsi(:, 1:1);
Bpsi_u = Bpsi(:, 2:end);
Dpsi_y = Dpsi(:, 1:1);
Dpsi_u = Dpsi(:, 2:end);

[n, ~] = size(A);
[ell, ~] = size(Apsi);
Ne = n+ell;

Ahat = [A, zeros(n,ell); Bpsi_y*C Apsi];
Bhat = [B; Bpsi_u+Bpsi_y*D];
Chat = [Dpsi_y*C, Cpsi];
Dhat = Dpsi_y*D+Dpsi_u;
%% Solve for rho with bisection using Yalmip
rhoub = 2.0;
rholb = 0.0;
rhotol = 10^(-4);

while (rhoub - rholb)>rhotol
    rho = (rhoub+rholb)/2;
    P = sdpvar(Ne);
    constraints = [P >= eye(Ne)]; %Want P PD, problem is scale invariant.
    H = sdpvar(N+1,N+1,'full');
    M = [zeros(N+1) H; H' zeros(N+1)];

    %constrains H to be doubly hyperdominant
    %off diagonal entries nonpositive
    for i=1:(N+1)
        for j=1:(N+1)
            if i~=j
                constraints = [constraints H(i,j)<=0];
            end
        end
    end

    %row and column sumns nonnegative
    for i=1:(N+1)
        constraints = [constraints sum(H(i,:))>=0];
        constraints = [constraints sum(H(:,i))>=0];
    end

    constraints = [constraints ([Ahat Bhat]'*P*[Ahat Bhat] - ...
              rho^2 *[eye(Ne) zeros(size(Bhat))]'*P*[eye(Ne) zeros(size(Bhat))]+ ...
              [Chat Dhat]'*M*[Chat Dhat])<=0];
    
    ops = sdpsettings('solver','mosek','verbose',0,'debug',0);
    diagnostics = optimize(constraints,0,ops);
    if diagnostics.problem == 0 %feasible
        rhoub = rho;
        P_star = value(P);
        M_star = value(M);
    else %infeasible
        rholb = rho;
    end
end
rho_star = rhoub;

% Z_star = [Ahat Bhat]'*P_star*[Ahat Bhat] - ...
%               rho^2 *[eye(Ne) zeros(size(Bhat))]'*P_star*[eye(Ne) zeros(size(Bhat))]+ ...
%               [Chat Dhat]'*M_star*[Chat Dhat];
% disp(eigs(Z_star))

    


