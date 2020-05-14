function [bub, R0, R1, R2] = nonasymptoticConvergenceRate(A,B,C,L,m,N)
% Variation of performance estimation problem to compute a bound b on the
% performance for the interconnection of the discxrete-time, linear system 
% system G:=(A(k),B(k),C(k),0) and function phi:
%   x(k+1) = A(k) x(k) + B(k) u(k)
%   y(k) = C(k) x(k)                     
%   u(k) = phi(y(k))
% Here phi is the gradient of a function f in S(m,L), i.e. f is strongly
% convex with parameter m and has gradient with Lipschitz bound L. 
% 
% Note that the inputs A, B, and C are passed as 3d arrays.
% the third dimension is time. 
%
% We let    eta = [x(1);u(1);...;u(N)]
%           x(0) = (R0 \kron I_d) eta
%           y(N) = (R1 \kron I_d) eta
%           z = Ly(1)-u(1);...;Ly(N)-u(1);u(1)-my(1);...;u(N)-my(N)]
%               = (R2 \kron I_d) eta 
% Define M_0 := b^2 R0'R0 - R1'R1
% 
% We find the smallest b such that at time N, || y(N+1) || <= b ||xg(1)||
% by solving(*) min: b 
%               subject to: M_0 - R2'[0 H; H' 0]R2>=0
% as provided in Theorem 1 in of:
%   1) Lee, Seiler. Finite Step Performance of First-order Methods Using Interpolation
%   Conditions Without Function Evaluations, arxiv. 2020.

%% Define linear maps from eta to x0,yN,z
Ng = size(A,1);

%R0 maps eta to x(0)
R0 = [eye(Ng) zeros(Ng,N)];

%R1 maps eta to y(N+1)
R1 = Prod(A); %term hitting x(1)
for i = 2:(N+1)
	R1 = [R1 Prod(A(:,:,i:end))*B(:,:,i-1)]; %term hitting u(i-1)
end
R1 = C(:,:,N+1) * R1; %Cg x(N+1) to y(N+1)

%R2 maps eta to z
R2 = [];
for i=1:N
    %ith iteration adds ith row to Rz which maps [xg(0); u(0);...;u(k-1)]
    %to y(i-1). Note that y() = Cg(:,:,i)*x(i). x(i) = sum_{ell=1:i-1}(
    %product_{k=ell+1, i-1}(A(:,:,k) Bg(:,:,ell) u_{})+ product_{k=0:i-1}(Ag(:,:,k))x(0)
    Rtmp = C(:,:,i)*Prod(A(:,:,1:i-1)); 
    for j = 1:(i-1)
        Rtmp = [Rtmp C(:,:,i)*Prod(A(:,:,(j+1):(i-1)))*B(:,:,j)];
    end
    Rtmp = [Rtmp zeros(1,N-i+1)];
    R2 = [R2; Rtmp];
end
R2 = [R2; zeros(N,Ng) eye(N)]; %pull the us off[xg(0); u(0);...;u(k-1)]
R2 = kron([L -1; -m 1],eye(N))*R2; %maps z to z_tilde

%% Solve (*) using yalmip by optimizing over b^2,H
b2 = sdpvar(1);
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

constraints = [constraints b2*(R0'*R0)-(R1'*R1)-(R2'*M*R2) >= 0];

ops = sdpsettings('solver','mosek','verbose',0,'debug',0);
optimize(constraints, b2, ops);
bub = sqrt(value(b2));