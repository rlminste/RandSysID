function z = hankelmult(x,c,r)
% Computes product of Hankel matrix with another matrix
%
% Inputs: 
%  x: matrix to be multiplied
%  c: last column of hankel matrix
%  r: first row of hankel matrix, reversed
%
% Output:
%  z: multiplication H*x, whether x is a matrix or vector
%
% - Written by Rachel Minster, 2020

[~,numcol] = size(x);

m = length(c);
z = zeros(m,numcol);

v = [c ; r(end:-1:2)']; 
fv = fft(v);

for j = 1:numcol
xc = x(:,j);
y = [xc(end:-1:1) ; zeros(m-1,1)]; 

zj = ifft(fft(y).*fv);
z(:,j) = zj(1:m);
end


end


