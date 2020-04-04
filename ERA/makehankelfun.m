function HF = makehankelfun(markov,s,l,m,tsignal)
%
% Sets up operator representing multiplication with block Hankel matrix
%
% Inputs: 
%   - markov: cell array of markov parameters
%   - s: number of block rows/columns of block Hankel matrix
%   - l: number of outputs
%   - m: number of inputs
%   - tsignal: 'transp' computes transpose multiplication,
%               'notransp' computes standard multiplication
%
% Outputs:
%   - HF: operator (funmat) representing multiplication with block Hankel
%   matrix

    N = 2*s-1; 
    
    % tsignal transposes the matrix implicitly.
    if nargin < 5, tsignal = 'notransp'; end

    % Define FV and FVT up here
    FV = zeros(N,l*m);
    FVT = zeros(N,l*m); 
    for i = 1:N
        h = reshape(markov{i},1,l*m); 
        ht = reshape(markov{i}',1,l*m); 
    
        j = mod(s+i-1,N)+1; 
        jt = mod(s+i-1,N)+1; 
    
        FV(j,:) = h;
        FVT(jt,:) = ht; 
    end 
    
    FV = fft(FV); 
    FVT = fft(FVT);
    
    if strcmp(tsignal,'notransp')
        HF = @(X,tflag) doHmult(X,tflag,s,s,l,m,FV,FVT);
    elseif strcmp(tsignal,'transp')
        HF = @(X,tflag) doHmult(X,tflag,s,m,l,FVT,FV); 
    else
        error("Transposition not indicated properly!"); 
    end
      
            
end

function HX = doHmult(X,tflag,s,l,m,FV,FVT)
% 
% Block Hankel multiplication
%
% Inputs: 
%   - X: matrix to be multiplied with
%   - tflag: 'transp' for transpose multiplication and
%           'notransp' for normal multiplication
%   - s: number of block rows/columns for block Hankel matrix
%   - l: number of outputs
%   - m: number of inputs
%   - FV: FFT of info from block Hankel matrix
%   - FVT: FFT of info from transpose block Hankel matrix
%
% Outputs:
%   - HX: product of block Hankel matrix and X

        nk = size(X,2); 
        N  = 2*s-1; 
        
        if strcmp(tflag,'transp')
            temp2 = l;
            l = m;
            m = temp2;
        elseif not(strcmp(tflag,'notransp'))
            error("Must indicate whether matrix is transposed!"); 
        end
        
        FX = zeros(s,m*nk);
        HX = zeros(s*l,nk);
        
            
        % FFT of the X data
        for idx = 0:(m-1)
            FX(s:-1:1,(idx*nk+1):((idx+1)*nk)) = X((idx+1):m:end,:);
        end
        FX = fft(FX,N); 

        % Compute the multiplication 
        for ixr = 0:(l-1)
            Z = zeros(N,nk); 
            for ixc = 0:(m-1)
                
                % Check which FFT set to multiply with
                if strcmp(tflag,'transp')
                    fv = FVT(:,ixc*l+ixr+1);
                else
                    fv = FV(:,ixc*l+ixr+1);
                end
                
                %Increment the fft-domain matvec product
                Z  = Z+FX(:,(ixc*nk+1):(ixc*nk+nk)).*fv; 
            end
            
            % Compute ifft and assign to output
            Z = ifft(Z); 
            HX((ixr+1):l:end,:) = Z(1:s,:);   
        end
end