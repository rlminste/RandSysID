function H = hankelize(markov,s,l,m)
% 
% Creates Hankel matrix from markov parameters

% Inputs:
%   markov: cell of markov parameters, each block l x m,  
%   s: number of block rows/columns
%   l: number of outputs
%   m: number of inputs
%
% - Written by Rachel Minster, 2020

H = zeros(s*l,s*m);

for i=1:s 
   inter = markov(i:s+i-1); 
   B = cell2mat(inter);
   H((i-1)*l+1:i*l,:) = B;
end

end
