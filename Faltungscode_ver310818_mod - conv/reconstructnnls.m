function [x, mu] = reconstructnnls(kernel,S,K,indices, y, convolve_flag ,opt)

%Full measurement matrix is NxMxL
N = size(kernel,1);
M = size(kernel,2);
L = length(indices);
class(kernel)

%Target Size
xdim = [N, M];

if strcmp(opt.method, 'lsqr') == 1
    %Reconstruction
    kernelF = fft2(kernel);
   % z = lsqr(@(X, flag) Afun(X, flag, kernelF,S,K, xdim, indices, convolve_flag), y, opt.tol, opt.max_iter);
  
      x = lsqnonneg(S,y);
    x = reshape(z(1:end-1), xdim);
    mu = z(end);
else
    error('no reconstruction method given');
end


function Y = Afun(X, flag, kernelF, S, K, xdim, indices, convolve_flag)
     if convolve_flag
         %X = X +0.05; 
%          for i= 1: size(X)
%              if X(i)>= 0
%                  X(i) = X(i);
%              else
%                  X(i)= 0;
%              end
%          end
         if strcmp(flag,'notransp')
            %Nontransposed 2DConvolution and Subsampling
            x = X(1:end-1);
            mu = X(end);
            %Y = conv2 (kernel, reshape(x, xdim));
            Y = ifft2(kernelF .* fft2(reshape(x, xdim)));
            Y = Y(indices) + mu;
            
        else
            %Transposed 2DConvolution and Subsampling
            tmp = zeros(xdim);
            tmp(indices) = X;
            %Y = conv2 (kernel, tmp);
            Y = (ifft2(conj(kernelF) .* fft2(tmp)));
            imagesc(Y)
            Y = [Y(:); sum(X(:))];
         end
     end