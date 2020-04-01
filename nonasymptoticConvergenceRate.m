function [bub, R0, R1, K] = nonasymptoticConvergenceRate(Ag,Bg,Cg,L,m,N,opt)
% Variation of performance estimation problem to compute a bound b on the
% performance for the interconnection of the discxrete-time, linear system 
% system G:=(Ag(k),Bg(k),Cg(k),0) and function phi:
%   xg(k+1) = Ag(k) xg(k) + Bg(k) u(k)
%   y(k) = Cg(k) xg(k)                     
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
%           z = [y(1);...;y(N);u(0);...;u(N)] = (R2 \kron I_d)eta 
%           z_tilde =
%               [Ly(1)-u(1);...;Ly(N)-u(1);u(1)-my(1);...;u(N)-my(N)]
%               = (K \kron I_d) eta 
% Define M_0 := b^2 R0'R0 - R1'R1
% 
% We find the smallest b such that at time N, || y(N+1) || <= b ||xg(1)||
% by solving(*) min: b 
%               subject to: M_0 - K'[0 H; H' 0]K>=0
% as provided in Theorem 2 in:
%   1) Lee, Seiler Finite Step Performance Guarantees for First Order Optimization
%   Algorithms: Revisited from a Control Theoretic Perspective, arxiv 2020
%
% Related Performance estimation problem is described in:
%   2) Taylor, Hendrickx, Glineur Smooth Strongly Convex Interpolation and 
%   Exact Worst-case Performance of First-order Methods, Mathematical
%   Programming 2017

%% Define linear maps from eta to x0,yN,z
Ng = size(Ag,1);

%R0 maps eta to xg(0)
R0 = [eye(Ng) zeros(Ng,N)];

%R1 maps eta to y(N+1)
R1 = Prod(Ag); %term hitting xg(1)
for i = 2:(N+1)
	R1 = [R1 Prod(Ag(:,:,i:end))*Bg(:,:,i-1)]; %term hitting u(i)
end
R1 = Cg(:,:,N+1) * R1; %maps xg(N+1) to y(N+1)

%R2 maps[xg(0); u(0);...;u(k-1)] to z_tilde
% util(N-1)]
R2 = [];
for i=1:N
    %ith iteration adds ith row to Rz which maps [xg(0); u(0);...;u(k-1)]
    %to y(i-1). Note that y() = Cg(:,:,i)*x(i). x(i) = sum_{ell=1:i-1}(
    %product_{k=ell+1, i-1}(A(:,:,k) Bg(:,:,ell) u_{})+ product_{k=0:i-1}(Ag(:,:,k))x(0)
    Rtmp = Cg(:,:,i)*Prod(Ag(:,:,1:i-1)); 
    for j = 1:(i-1)
        Rtmp = [Rtmp Cg(:,:,i)*Prod(Ag(:,:,(j+1):(i-1)))*Bg(:,:,j)];
    end
    Rtmp = [Rtmp zeros(1,N-i+1)];
    R2 = [R2; Rtmp];
end
R2 = [R2; zeros(N,Ng) eye(N)]; %pull the us off[xg(0); u(0);...;u(k-1)]
K = kron([L -1; -m 1],eye(N))*R2; %maps z to z_tilde

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

constraints = [constraints b2*(R0'*R0)-(R1'*R1)-(K'*M*K) >= 0];

ops = sdpsettings('solver','mosek','verbose',0,'debug',0);
optimize(constraints, b2, ops);
bub = sqrt(value(b2));