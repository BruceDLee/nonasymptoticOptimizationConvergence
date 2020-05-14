%This file tests code to find a worst case nonasymptotic trajectory
%by defining an algorithm and function class SmL, and checking
%that the generated trajectory is interpolable in SmL. It is then verified 
%that it attains the worst case bound.
L = 10;
m = 1;
N = 12; 

%define time varying algorithm
alpha = zeros(N);
beta = zeros(N);

%let algorithm be time invariant heavy ball optimized for quadratics
for i=1:N
    alpha(i) =  4/(sqrt(L)+sqrt(m))^2; 
    beta(i) =  ((sqrt(L)-1)/(sqrt(L)+1))^2;
end
    
Ag = zeros(2,2,N);
Bg = zeros(2,1,N);
Cg = zeros(1,2,N+1);

for i=1:N
    Ag(:,:,i) = [1+beta(i), -beta(i); 1, 0];
    Bg(:,:,i) = [-alpha(i); 0];
    Cg(:,:,i) = [1 0];
end

Cg(:,:,N+1) =  [1 0];

[bub, blb, B_star, R2] = nonasymptoticWorstCaseTraj(Ag,Bg,Cg,L,m,N);
%B_star is [x(0)' (reshaped to R^(d \times n);u(0)';u(1)';...;u(N-1)']
%R2 maps eta to z

%compute full trajectory from x0 and inputs
[~,d] = size(B_star);
xs = zeros(N+1,2*d);
xs(1,:) = [B_star(1,:), B_star(2,:)];
us = B_star(3:end,:);

ys = zeros(N+1, d);
ys(1,:) = (kron(Cg(:,:,1), eye(d))*xs(1,:)')';

for i = 1:N
    xs(i+1,:) = (kron(Ag(:,:,i),eye(d))*xs(i,:)'+kron(Bg(:,:,i),eye(d))*us(i,:)')';
    ys(i+1,:) = (kron(Cg(:,:,i+1), eye(d))*xs(i+1,:)')';
end
%%
bub
blb
disp('||y_N||/||x_0||:');
disp(norm(ys(end,:))/norm(xs(1,:)));

%% Used for debugging
% The code segment below checks that the trajectory is feasible
% H = zeros(N+1,N+1,factorial(N+1));
% Ps = perms(1:N+1);
% I = eye(N+1);
% J = [zeros(1, N); eye(N)];
% 
% for i = 1:(factorial(N+1))
%     H(:,:,i) = I - I(Ps(i,:),:);
% end
% 
% %zs = (kron([L, -1; -m, 1], eye(N*d)))*[reshape(ys(1:end-1,:), [N*d 1]); reshape(us,[N*d 1])];
% zs = kron(Rz, eye(d))*[reshape(xs(1, :)', [2*d 1]); reshape(us', [N*d 1])];
% 
% %check trajectory
% disp('Check trajectory:')
% lb = Inf;
% for i=1:factorial(N+1)
%     lb = min(zs'*kron([zeros(N), J'*H(:,:,i)*J; (J'*H(:,:,i)*J)',zeros(N)], eye(d))*zs, lb);
% end
% disp('Lower bound on [ytil][0; I]^T(I - Pk)[0; I]^T[util] (nonnegative => interpolable):');
% disp(lb);

