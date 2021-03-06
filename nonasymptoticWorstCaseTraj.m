function [bub, blb, B_star, R2] = nonasymptoticWorstCaseTraj(A,B,C,L,m,N)
%% Description 
% again we have
%G:=(A(k),B(k),C(k),0) and function phi:
%   x(k+1) = A(k) x(k) + B(k) u(k)
%   y(k) = C(k) x(k)                     
%   u(k) = phi(y(k))
% Here phi is the gradient of a function f in S(m,L), i.e. f is strongly
% convex with parameter m and has gradient with Lipschitz bound L. 
% 
% Note that the inputs Ag, Bg, and Cg are passed as 3d tensors.
% the third dimension is time. 
%
% We let    eta = [xg(1);u(1);...;u(N)]
%           xg(0) = (b R0 \kron I_d) eta
%           y(N) = (R1 \kron I_d) eta
%           z = [Ly(1)-u(1);...;Ly(N)-u(1);u(1)-my(1);...;u(N)-my(N)]
%               = (R2 \kron I_d) eta 
% Define M_0 := b^2 R0'R0 - R1'R1
% 
% We find the smallest b such that at time N, || y(N+1) || <= b ||xg(1)||
% by solving(*) b_* = min: b 
%               subject to: M_0 - K'[0 H; H' 0]K>=0
% as provided in Theorem 2 in:
%   1) Finite Step Performance of First-order Methods Using Interpolation
%   Conditions Without Function Evaluations. arxiv. 2020
%
% We then find a sequence attaining the nonasymptotic worst case behavior.
%
% We use the fact that for any b<b^*, the lmi,
% R1'R1 - b^2*R0'R0-R2[0 H; H' 0]R2<= 0 is infeasible for all 
% doubly hyperdominant H. 
% Then the optimal solution to the above problem
%to the problem
%   max_{nu, H doubly hyperdominant}: nu
%   subject to: R1'R1 - b^2R0'R0 - R2'[0 H; H' 0] R2 <= nuI
%is negative. 
%
% As dicussed in 1), a solution H_* may be decomposed into a conic
% combinations of cyclic monotonicity constraint matrices,
% which provides a solution (nu_star, lambda_star) to 
%   min_{nu, lambda}: nu
%   subject to: R1'R1 - b^2R0'R0 - R2' (sum lambda_i J'(I-P-i)J) R2  <= nuI
% where J = [0; I];
%
% Strong duality holds for this problem, so the dual:
%min_{G}: trace((R1'R1 - b^2 R0'R0)G)
%subject to: G >= 0
%            trace(G) == 1
%            trace(G Rz'MiRz) >=0 for i with lambda_{*,i}>0
%also has an optimal objective less than 0. As dicussed in 1),
% if G_* is unique, it can be factored to obtain a feasible
% trajectory, R = [x(0); u(0);...;u(N-1)], which violates ||yN|| <= b||x0||
%
% The matrix B_star that is returned contains an x0 and set of N u values
% that create an interpolable trajectory, and serve as a lower bound on
% performance.
%%
%get upper and lower bounds on convergence rate as well as helper matrices
[bub, R0,R1,R2] = nonasymptoticConvergenceRate(A,B,C,L,m,N);
Ng = size(A,1);
blb = bub-10^(-3);
M0 = blb^2*R0'*R0 - R1'*R1; 
[k,l] = size(M0);
%% solve dual problem
[~, H_star, Z_star] = solvedual(M0, R2);
%% solve primal problem 
G = sdpvar(k,l);
constraints = [G >= 0, trace(G) == 1];
I = eye(N+1);
%% (all constraints) (necessary if (nu_*, lambda_*) degenerate)
% Ps = perms(1:N+1);
% P = zeros(N+1,N+1, size(Ps,1));
% for i =1:size(Ps,1)
%    P(:,:,i) = I(Ps(i,:),:);
% end
% disp(size(P))
%% (Birkhoff constraints only)
%H_star is dhd matrix solution to problem equivalent to dual
%make H_hat dhd zero excess
H_hat = [-sum(H_star, 1); H_star];
H_hat = [-sum(H_hat,2),H_hat]; 

%Permutation matrices corresponding to lambda_{*,i}>0
[P, ~] = dhdecomp(H_hat); 
%% check nondegeneracy of (nu_*, lambda_*)
[V,D] = eig(Z_star);
s = length(diag(D(D>10^(-6))));
Q1 = V(s+1:end,:);
n_s = size(Q1, 1);

Bases = zeros(n_s,n_s,size(P,3)+1);
Bases(:,:,1) = eye(n_s);

Bases_vectorized = zeros(n_s*(n_s+1)/2);
function v = uppertrianglevectorize(B)
    if size(B,1) < 2
        v = B;
    else 
        v = [B(:,1); uppertrianglevectorize(B(2:end,2:end))];
    end
end
Bases_vectorized(:,1) = uppertrianglevectorize(Bases(:,:,1));

H = zeros(N+1,N+1, size(P,3));
J = [zeros(1, N); eye(N)];
for i = 1:(size(P,3))
    H(:,:,i) = I - P(:,:,i);
    Bases(:,:,i+1) = Q1*R2'*[zeros(N),J'*H(:,:,i)*J; ...
        (J'*H(:,:,i)*J)',zeros(N)]*R2*Q1';
    Bases_vectorized(:,i+1) = uppertrianglevectorize(Bases(:,:,i+1));
end

[~,S,~] = svd(Bases_vectorized);

if n_s*(n_s+1)/2~=length(S(S>10^(-6)))
    disp('Nondegeneracy conditions not satisfied, solution invalid. Use all constraints') 
end

%% 
for i = 1:(size(P,3))
    Q(:,:,i) = I - P(:,:,i);
    constraints = [constraints, trace(G*R2'*[zeros(N),J'*Q(:,:,i)*J; ...
        (J'*Q(:,:,i)*J)',zeros(N)]*R2) >= 0];
end
ops = sdpsettings('solver','mosek','verbose',0,'debug',0);
optimize(constraints, trace(M0*G), ops);
G_star = value(G);
 
%%
%Now we rank factorize G_star
[u, s, ~] = svd(G_star);
diagS = diag(s);
nonZeroSVs = diagS(diagS >10^(-6));
d = length(nonZeroSVs);
B_star = u*sqrt(s);
B_star = B_star(:, 1:d);

end