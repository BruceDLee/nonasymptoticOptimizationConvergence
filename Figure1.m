% Generates Figure 1 in 
%   1) Lee, Seiler Finite Step Performance Guarantees for First Order Optimization
%   Algorithms: Revisited from a Control Theoretic Perspective, arxiv 2020

m = 1;
Ls = [10, 30];
for L=Ls
    Ndat=4:12;
    JJ=zeros(size(Ndat));
    for ii=1:numel(Ndat)
        N=Ndat(ii);
        
        %Define N steps of heavy ball algorithm optimized for
        %quadratics(m,L):
        alpha = zeros(N);
        beta = zeros(N);

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

        %compute upper bound on N-step performance bound
        [b,~,~,~] = nonasymptoticConvergenceRate(Ag,Bg,Cg,L,m,N);

        JJ(ii)=b;
    end

    rhos = zeros(3,1);
    %compute asymptotic bound defined in 1)
    for T=1:3
        cyclic_pointwise_IQC_Yalmip
        rhos(T) = rho_star;
    end

    close all;
    plot(Ndat,JJ,'b-*',Ndat,JJ(1)*rhos(1).^(Ndat-Ndat(1)),'k', ...
        Ndat,JJ(1)*rhos(2).^(Ndat-Ndat(1)),'g',Ndat,JJ(1)*rhos(3).^(Ndat-Ndat(1)),'r')
    legend('Finite', 'T=1', 'T=2', 'T=3')
    set(gca,'FontSize',30)
    xlabel('N')
    ylabel('b_*')
    %garyfyFigure
    %print(gcf, "Figure1_L"+num2str(L)+".pdf", '-dpdf', '-fillpage')
    %p = get(gcf,'Position');
    %set(gcf,'Position',[p(1) p(2) p(3) p(4)*3/5]);
    %myprint(gcf,"Figure1_L"+num2str(L))
end
   
