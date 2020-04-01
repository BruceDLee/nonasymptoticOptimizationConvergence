%compare the bound produced from performance estimation problem package:
%   1) Taylor, Hendrickx, Glineur Smooth Strongly Convex Interpolation and 
%   Exact Worst-case Performance of First-order Methods, Mathematical
%   Programming 2017
%to that found in
%   2) Lee, Seiler Finite Step Performance Guarantees for First Order Optimization
%   Algorithms: Revisited from a Control Theoretic Perspective, arxiv 2020

% consider worst case N-step performance of the algorithm with decaying
% stepsize and momentum parameter on SmL using the approaches in each
% reference. The displayed outputs are the N-step performance estimation
% bound defined in 2) squared.

L = 10;
m = 1;

P = pep();
param.mu = m;
param.L = L;

F=P.DeclareFunction('SmoothStronglyConvex',param);

yminus = P.StartingPoint();
y0=P.StartingPoint();
[ys,fs]=F.OptimalPoint();
P.InitialCondition((y0-ys)'*(y0-ys)<=0); %+(yminus-ys)'*(yminus-ys)<=1);
%sets \|x_0\| \leq 1

N=15;
alpha = zeros(N);
beta = zeros(N);

for i=1:N
    alpha(i) = 1/param.L/i; %decaying stepsize
    beta(i) = 0.01;
end

y_ = yminus;
y = y0;
for i=1:N
    ytemp = y;
    y =(1+beta(i))*y-beta(i)*y_-alpha(i)*F.gradient((1+beta(i))*y-beta(i)*y_);
    y_ = ytemp;
end
yN_ = y_;
yN = y;

% Set up the performance measure
P.PerformanceMetric((yN-ys)'*(yN-ys));

% Solve the PEP and evaluate output
P.solve()
disp('1 provides:')
disp(double(yN-ys)'*double(yN-ys));

%evaluate approach from 2)
Ag = zeros(2,2,N);
Bg = zeros(2,1,N);
Cg = zeros(1,2,N+1);

for i=1:N
Ag(:,:,i) = [1+beta(i), -beta(i); 1, 0];
Bg(:,:,i) = [-alpha(i); 0];
Cg(:,:,i) = [1+beta(i), -beta(i)];
end

Cg(:,:,N+1) = [1+beta(N), -beta(N)];

[bub, blb,h,R1, R0, Rz] = nonasymptoticConvergenceRate(Ag,Bg,Cg,L,m,N);
disp('2 provides:')
disp(bub^2)