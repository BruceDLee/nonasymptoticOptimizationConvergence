%% compare asymptotic results
%
% A demonstration that Theorem 2 from 
%   Lee, Seiler.Finite Step Performance of First-order Methods Using Interpolation
%   Conditions Without Function Evaluations. arxiv, 2020.
% attains tighter convergence rate bounds on heavy ball algorithm than
% weighted of by one IQC from
%   Lessard, Recht, Packard. Analysis and Design of Optimization Algorithms 
%   via Integral Quadratic Constraints. SIAM, 2016.
%   
m = 1;
Ls = 1:50;
Ts = [1, 2, 3];
rhos = zeros(size(Ts,2)+1,size(Ls,2));

for enn = 1:size(Ls,2)
    for emm = 1:size(Ts,2)
        T = Ts(emm);
        L = Ls(enn);
        cyclic_pointwise_IQC
        rhos(emm,enn) = rho_star;
    end
    alpha = 4/(sqrt(L)+sqrt(m))^2;
    beta = ((sqrt(L)-sqrt(m))/(sqrt(L)+sqrt(m)))^2;
    A = [1+beta -beta; 1 0];
    B = [-alpha;0];
    C = [1+beta -beta];
    D = 0;
    G = ss(A,B,C,D,1);
    rhos(4, enn) = ratebound(G, L,m,1);
end

%% comparison plot
for enn = 1:size(Ls,2)
    L = Ls(enn);
    alpha = 4/(sqrt(L)+sqrt(m))^2;
    beta = ((sqrt(L)-sqrt(m))/(sqrt(L)+sqrt(m)))^2;
    A = [1+beta -beta; 1 0];
    B = [-alpha;0];
    C = [1, 0];
    D = 0;
    G = ss(A,B,C,D,1);
    rhos(4, enn) = ratebound(G, L,m,2);
end
plot(Ls, rhos(1,:), 'b', Ls, rhos(2,:), 'c', Ls, rhos(3,:), 'k', Ls, rhos(4,:), 'r--')
legend('T=1', 'T=2', 'T=3', 'WgtedOffbyOne')
xlabel('L/m')
ylabel('\rho')