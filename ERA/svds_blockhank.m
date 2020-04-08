function y = svds_blockhank(x,tflag,T,l,m,s)
%
% fits svds function with block hankel multiplication
%
% Inputs: 
%       x: vector to be multiplied
%       tflag: indicates whether multiplying H or H'
%           - 'notransp' multiplies with H
%           - 'transp' multiplies with H'
%       T: matrix representing block hankel matrix H
%       l: number of outputs
%       m: number of inputs
%       s: number of block rows/columns of H
%
% Outputs: y = Hx or H'x
%
% - Written by Arvind K. Saibaba and Rachel Minster, 2020

    if strcmp(tflag,'notransp')
        y = blockhankmult(T,x,s,l,m,0);
    elseif strcmp(tflag,'transp')
        y = blockhankmult(T,x,s,l,m,'transpose');
    end
    
end