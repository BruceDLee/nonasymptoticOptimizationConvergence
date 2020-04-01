function [P,a,kmax] = birkhoff(M)
% birkhoff  Birkhoff-von Neumman decomposition for a doubly stochastic matrix
%
% [P,alpha]=birkhoff(M) returns a Birkhoff-von Neumman decomposition for 
% an n-by-n doubly stochastic matrix M. This returns P as an n-by-n-by-k
% array and alpha as a k-by-1 vector such that:
%    M = alpha(1)*P(:,:,1) + ... + alpha(k)*P(:,:,k)
% where each P(:,:,i) is a permutation matrix. In addition the entries of
% alpha are non-negative and sum to 1. 
%
% % Example: Create random doubly stochastic matrix 
% n=5;  I=eye(n); M = zeros(n); 
% for i=1:100, 
%    p = randperm(n);
%    M=M+abs(randn)*I(:,p); 
% end; 
% M = M/sum(M(:,1));
%
% % Compute decomposition and verify result
% [P,alpha] = birkhoff(M);  
% M0=zeros(n); 
% for i=1:numel(alpha), 
%    M0=M0+alpha(i)*P(:,:,i); 
% end
% [norm(M0-M) abs(1-sum(alpha))]
%
% References
% [1] Birkhoff, Tres observacionea sobre el algebra lineal, Univ. Nac. 
%  Tucuman Rev, Ser. A, Vol. 5, p. 147-150, 1946.
% [2] von Neumann, A certain zero-sum two-person game equivalent to the 
%  optimal assignment problem, Contributions to the Theory of Games, 
%  Princeton University Press, Vol. 2, pp. 5-12, 1953.
% [3] Marshall and Olkin, Inequalities theory of majorization and its 

% Tolerance for detecting zero entries
tol = 100*eps;

% Verify that M is a square, real matrix
if ~( ismatrix(M) && size(M,1)==size(M,2) && isreal(M) )
    error('M must be a square, real matrix')
end
n = size(M,1);

% Verify that M is doubly stochastic
if ~all(M(:)>=0) 
    error('All entries of M must be non-negative.');
end
if ~( max(abs(sum(M,1)-1)<tol) && max(abs(sum(M,2)-1)<tol) )
    error('All rows and columns of M must sum to 1.')
end

% Birkhoff Algorithm
% DMPERM below uses the Dulmage and Mendelsohn method
I = eye(n);
kmax = n^2-2*n+2;

P = zeros(n,n,kmax);
a = zeros(kmax,1);
k = 0;
while max(abs(M(:)))>tol
    k = k+1;
    p = dmperm(M);
    Pk = I(:,p);
    ak= min( M(Pk~=0) ); 
    M=M-ak*Pk;
    
    M( abs(M)<tol ) = 0;
    P(:,:,k) = Pk;
    a(k) = ak;
end

if k>kmax
    warning('Iteration took more iterations than theoretical upper bound.')
end
P = P(:,:,1:k);
a = a(1:k);

%% DETAILS
% A Birkhoff von-Neumann decomposition [1,2,3] can be done with
% k<= kmax:=n^2-2*n+2 by the Marcus-Ree theorem [3,4]. This is the
% tightest bound that holds over all doubly stochastic matrices [3,5].
% Finding a decomposition with the smallest k is strongly NP-complete [6].
% This function uses the Birkhoff algorithm to construct a decomposition.
% This algorithm finds a decomposition with no more than kmax terms
% [7,8]. A key step of this algorithm requires the calculation of a 
% maximum matching for a bipartite graph. This step is performed using 
% the Dulmage-Mendelsohn method as implemented in DMPERM.
%
% References
% [1] Birkhoff, Tres observacionea sobre el algebra lineal, Univ. Nac. 
%  Tucuman Rev, Ser. A, Vol. 5, p. 147-150, 1946.
% [2] von Neumann, A certain zero-sum two-person game equivalent to the 
%  optimal assignment problem, Contributions to the Theory of Games, 
%  Princeton University Press, Vol. 2, pp. 5-12, 1953.
% [3] Marshall and Olkin, Inequalities theory of majorization and its 
%  applications, Academic Press, 1979.
% [4] Marcus and Ree, Diagonals of doubly stochastic matrices, The 
%  Quarterly Journal of Mathematics, Vol. 10, No. 1, p. 296–302, 1959.
% [5] Farahat and Mirsky, Permutation endomorphisms and refinement of a 
%  theorem of Birkhoff, Mathematical Proceedings of the Cambridge 
%  Philosophical Society, vol. 56, no. 4, p.322-328, 1960.
% [6] Dufossé and Uçar, Notes on Birkhoff–von Neumann decomposition of
%  doubly stochastic matrices, Linear Algebra and its Applications, 
%  vol. 497, p. 108–115, 2016.
% [7] Johnson, Dulmage, and Mendelsohn, On an algorithm of G. Birkhoff
%  concerning doubly stochastic matrices, Canadian Math. Bulletin, vol. 3 
%  p. 237-242, 1960.
% [8] Brualdi, Notes on the Birkhoff algorithm for doubly stochastic 
%  matrices, Canadian Math. Bulletin, vol. 25, no. 2, p. 191-199, 1982.
% [9] Dulmage and Mendelsohn, Coverings of bipartite graphs, Canadian 
%  Journal of Mathematics, vol. 10, p.517–534, 1958.
