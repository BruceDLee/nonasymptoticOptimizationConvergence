function [P,b] = dhdecomp(M)
% dhdecomp  Decomposition for a doubly hyperdominant matrix with no excess
%
% [P,beta]=dhdecomp(M) returns decomposes an n-by-n doubly hyperdominant 
% matrix M with no excess. This returns P as an n-by-n-by-k array and 
% beta as a k-by-1 vector such that:
%    M = beta(1)*(I-P(:,:,1)) + ... + beta(k)*(I-P(:,:,k))
% where each P(:,:,i) is a permutation matrix. In addition the entries of
% beta are non-negative. 
%
% % Example: Create random doubly hyperdominant matrix with no excess
% n=5;  I=eye(n); S = zeros(n); 
% for i=1:100, 
%    p = randperm(n);
%    S=S+abs(randn)*I(:,p); 
% end; 
% S = S/sum(S(:,1));
% r = abs(randn);
% M = r*(I-S);
%
% % Compute decomposition and verify result
% [P,beta] = dhdecomp(M);  
% M0=zeros(n); 
% for i=1:numel(beta), 
%    M0=M0+beta(i)*(I-P(:,:,i)); 
% end
% norm(M0-M)
%
%
% References
% [1] Willems, The Analysis of Feedback Systems, MIT Press, 1971.
% (See Section 3.5 and the proof of Theorem 3.7)

% Tolerance for detecting zero entries
tol = 100000*eps;

% Verify that M is a square, real matrix
if ~( ismatrix(M) && size(M,1)==size(M,2) && isreal(M) )
    error('M must be a square, real matrix')
end
n = size(M,1);

% Verify that M is doubly hyperdominant
Mdiag = diag(M);
if ~all(Mdiag>=0)
    error('Diagonal entries of M must be non-negative.')
end
if ~all( M-diag(Mdiag)<=0 )
    error('Off-diagonal entries of M must be non-positive.')
end
if ~( all(abs(sum(M,1))<=tol) && all(abs(sum(M,1))<=tol) )
    error('The sum along each row and column of M must be zero.')
end

% Decompose M = r*[I-S] where S:=(r*I-M)/r is doubly stochastic
% (As in the proof of Theorem 3.7 of Willems)
r = max( abs(M(:)) );
S = eye(n)-M/r;

% Compute a Birkhoff-von Neumman Decomposition:
%    S = a(1)*P(:,:,1) + ... + a(k)*P(:,:,k)
% Each P(:,:,i) is a perm. matrix, sum(a)=1, and each a(i)>=0.
[P,a]=birkhoff(S);

% Convert to a conic combination:
%   M = b(1)*[I-P(:,:,1)] + ... + b(k)*[I-P(:,:,k)]
% where b(i) = r*a(i) and each b(i)>=0.
b = r*a;

%% COMMENT
% If M is doubly-hyperdominant with excess then an additional row and 
% column can be appended to create an (n+1)-by-(n+1) matrix Mtil with
% no excess:
%    % Note: sum(rowsum)=sum(colsum) and either can be used for Mtil(1,1)
%    colsum = sum(M,1);
%    rowsum = sum(M,2);
%    Mtil = [sum(rowsum) -colsum; rowsum M];
