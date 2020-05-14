function prod = Prod(A)

n = size(A,1);
m = size(A,3);
if m == 0
    prod = eye(n);
else
    prod = Prod(A(:,:,2:end))*A(:,:,1);
end