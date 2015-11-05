function getEigenValues_FD4

stencil = (1/12)*[-1 16 -30 16 -1];
n=32;   % size of Matrix
x0 = zeros(1,n);
x0(1:5) = stencil;

y = circshift(x0,[0, -2]);

M = zeros(n,n);
for i = 1:n
   M(i,:) =y;
   y = circshift(y,[0, 1]);
end

evs = eig(M);
d = eigs(M);