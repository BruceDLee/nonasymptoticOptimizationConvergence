%compare the performance estimation problem package provided by 
%   Taylor, Hendrickx, Glineur Smooth Strongly Convex Interpolation and 
%   Exact Worst-case Performance of First-order Methods, Mathematical
%   Programming 2017
%to our method.
%their code is available here at:
%   https://github.com/AdrienTaylor/Performance-Estimation-Toolbox
%and is required. Test on heavy ball optimized for quadratics.
% the time to compute the bound and worst case trajectory is also computed.

L = 10;
m = 1;

P = pep();
param.mu = m;
param.L = L;

F=P.DeclareFunction('SmoothStronglyConvex',param);
Ns = 2:10;

times = zeros(2, size(Ns, 2));
bounds = zeros(2, size(Ns,2));

yminus = P.StartingPoint();
y0=P.StartingPoint();
[ys,fs]=F.OptimalPoint();
P.InitialCondition((y0-ys)'*(y0-ys)+(yminus-ys)'*(yminus-ys)<=1);

for j = 1:size(Ns,2)
    
    N = Ns(j);
    alpha = zeros(N+1);
    beta = zeros(N+1);

    for i=1:N
        alpha(i) = 4/(sqrt(L)+sqrt(m))^2; 
        beta(i) = ((sqrt(L)-sqrt(m))/(sqrt(L)+sqrt(m)))^2;
    end

    y_ = yminus;
    y = y0;
    for i=1:N
        ytemp = y;
        y =(1+beta(i))*y-beta(i)*y_-alpha(i)*F.gradient(y);
        y_ = ytemp;
    end
    yN_ = y_;
    yN = y;

    % Set up the performance measure
    P.PerformanceMetric((yN-ys)'*(yN-ys));

    % Solve the PEP and evaluate output
    tic
    P.solve()
    times(1,j) = toc;
    bounds(1,j) = sqrt(double(yN-ys)'*double(yN-ys));

    %evaluate our approach
    Ag = zeros(2,2,N);
    Bg = zeros(2,1,N);
    Cg = zeros(1,2,N+1);

    for i=1:N
        Ag(:,:,i) = [1+beta(i), -beta(i); 1, 0];
        Bg(:,:,i) = [-alpha(i); 0];
        Cg(:,:,i) = [1, 0];
    end

    Cg(:,:,N+1) = [1, 0];
    tic
    [bub, blb, B_star, K] = nonasymptoticWorstCaseTraj(Ag,Bg,Cg,L,m,N);
    [~,d] = size(B_star);
    xs = zeros(N+1,2*d);
    xs(1,:) = [B_star(1,:), B_star(2,:)];
    us = B_star(3:end,:);

    zs = zeros(N+1, d);
    zs(1,:) = (kron(Cg(:,:,1), eye(d))*xs(1,:)')';

    for i = 1:N
        xs(i+1,:) = (kron(Ag(:,:,i),eye(d))*xs(i,:)'+kron(Bg(:,:,i),eye(d))*us(i,:)')';
        zs(i+1,:) = (kron(Cg(:,:,i+1), eye(d))*xs(i+1,:)')';
    end
    times(2,j) = toc;
    bounds(2,j) = bub;
end

%%
close all
plot(Ns, bounds(1,:), 'r', 'LineWidth', 4)
hold on
plot(Ns, bounds(2,:), 'k--', 'LineWidth', 4)
xlabel('N')
ylabel('bound')
legend('PEPtoolbox','No function evals approach')
figure()
plot(log(Ns), log(times(1,:)), 'r', 'LineWidth', 4)
hold on
plot(log(Ns), log(times(2,:)), 'k--', 'LineWidth', 4)
xlabel('log(N)')
ylabel('log(time)')
legend('PEPtoolbox', 'No function evals approach')