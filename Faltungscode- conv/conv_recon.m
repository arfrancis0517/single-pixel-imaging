function x = conv_recon(kernel, indices, y)

%Full measurement matrix is NxMxL
N = size(kernel,1);
M = size(kernel,2);
L = length(indices);

%Target Size
xdim = [N, M];

%FlexBox - new Varient
main = flexBox;
main.params.maxIt = 50;
main.params.relativePathToMEX = '../flexBox_CPP/source/build/bin';
main.params.tryCPP = 0;
main.params.checkError = 1;
main.params.showPrimals = 0;

op = concatOperator(...
		subsamplingOperator(indices, xdim), ...
		convolutionOperator(kernel, xdim),...
		'composition');
	
nX = main.addPrimalVar(xdim);

main.addTerm(L2dataTermOperator(1/L, op, y), nX);
main.addTerm(L2gradient(0.01, xdim), nX);
main.addTerm(nonNegativityConstraint(xdim), nX);

for i = 1:1000
    main.runAlgorithm();
    x = main.getPrimal(nX);
    imagesc(x);
    drawnow;
end

end

%Reconstruction MATLAB Variation

% %Function Handle for Matlab
% function Y = Afun(X, flag)
%     if strcmp(flag,'notransp')
%         %Nontransposed 2DConvolution and Subsampling
%         Y = ifft2(kernelF .* fft2(reshape(X, xdim)));
%         Y = Y(indices);
%     else
%         %Transposed 2DConvolution and Subsampling
%         tmp = zeros(xdim);
%         tmp(indices) = X;
%         Y = vec(ifft2(conj(kernelF) .* fft2(tmp)));
%     end
% end
% 
% %Reconstruction
% x = lsqr(@Afun, y);
% x = reshape(x, xdim);
        
%Old Flexbox
% function Y = AfunNotransp(X)
% 	Y = ifft2(kernelF .* fft2(reshape(X,xdim)));
% 	Y = Y(indices);
% end
% 
% function Y = AfunTransp(X)
% 	tmp = zeros(xdim);
%  	tmp(indices) = X;
% 	Y = vec(ifft2(conj(kernelF) .* fft2(tmp)));
% end
% 
% op = functionHandleOperator(@AfunNotransp, @AfunTransp, xdim, norm(vec(fft2(kernel)),Inf));
        

