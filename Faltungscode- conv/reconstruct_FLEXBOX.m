function [x, mu] = reconstruct(kernel,S,K,indices, y, convolve_flag ,opt)

%Full measurement matrix is NxMxL
N = size(kernel,1);
M = size(kernel,2);
L = length(indices);
class(kernel)

%Target Size
xdim = [N, M];



if strcmp(opt.method, 'flexbox') == 1
    %FlexBox
    main = flexBox;
    main.params.relativePathToMEX = '\flexBox_CPP\build\bin';
    main.params.tryCPP = 1;
    main.params.maxIt = opt.flexbox_internal_iter;
    main.params.checkError = 1;
    main.params.showPrimals = 0;
    main.params.verbose = 1;
    op = concatOperator(...
            subsamplingOperator(indices, xdim), ...
            convolutionOperator(kernel, xdim),...
            'composition');

    nX = main.addPrimalVar(xdim);
    nMu = main.addPrimalVar(1);
    
 for it = 1:opt.max_iter   
    if 1
        kernelF = fft2(kernel);
       % z = lsqr(@(X, flag) Afun(X, flag, kernelF,S,K, xdim, indices, convolve_flag), y, opt.tol, opt.max_iter);
         z = lsqr(S,y,opt.tol,opt.max_iter);
        main.x{nX} = z(1:end-1);
        main.x{nMu} = z(end);
    end
    
    if opt.flexbox_l1grad > 0
        main.addTerm(L1gradientIso(opt.flexbox_l1grad, xdim), nX);
    end
    if opt.flexbox_l2grad > 0
        main.addTerm(L2gradient(opt.flexbox_l2grad, xdim), nX);
    end
    if opt.flexbox_nnls == 1
        main.addTerm(nonNegativityConstraint(xdim), nX);
    end
    if isfield(opt, 'flexbox_wname') && strcmp(opt.flexbox_wname, 'none') == 0
        wname = opt.flexbox_wname;

        % measure the maximal depth we can archieve with this wavelet
        % transform
        lvl = wmaxlev(xdim, wname);
                
        % assume structural component ist always the same, also get the
        % output size oif the wavelet transform
        [C, S] = wavedec2(zeros(xdim), lvl, wname);
        wdim = [numel(C), 1];
        clear C;
        
        %Estimate the norm of the Wavelet operator
        [~, opNorm, ~] = svds(@(X, flag) Wfun(X, flag, lvl, S, wname, xdim), [wdim(1), N*M], 1, 'largest');
        
        %Give over a function handle for the wavelet transform
        waveop = functionHandleOperator(@(x) WaveletHandle(x, lvl, wname, xdim), @(x) WaveletHandleT(x, S, wname, xdim), xdim, opNorm);
        
        %Finally regularize on the L1 norm of Wx
        %main.addTerm(L1dataTermOperator(opt.flexbox_wfactor, waveop, zeros(wdim)), nX); 
    end
    %main.addTerm(L2dataTermOperator(1/L, {op, ones(L,1)}, y), [nX, nMu]);
        main.runAlgorithm();
        x = main.getPrimal(nX);
        mu = main.getPrimal(nMu);
%     tstart = tic();

%     for it = 1:opt.max_iter
%         main.runAlgorithm();
%         % x= zeros(3072:1);
%         x = main.getPrimal(nX);
%         mu = main.getPrimal(nMu);
%        
%         if isfield(opt, 'callback')
%             opt.callback(x, mu, it, 0);
%         end
% 
%         time = toc(tstart);
%         if time > opt.max_time
%             break;
%         end
 end
%     
elseif strcmp(opt.method, 'lsqr') == 1
    %Reconstruction
    kernelF = fft2(kernel);
    z = lsqr(@(X, flag) Afun(X, flag, kernelF,S,K, xdim, indices, convolve_flag), y, opt.tol, opt.max_iter);
    x = reshape(z(1:end-1), xdim);
    mu = z(end);
else
    error('no reconstruction method given');
end

end

function Y  = Wfun(X, flag, lvl, S, wname, xdim)
    if strcmp(flag, 'notransp')
        %Nontransposed wavelet transform
        Y = WaveletHandle(X, lvl, wname, xdim);
        Y = Y(:);
    else
        Y  = WaveletHandleT(X, S, wname, xdim);
        Y = Y(:);
    end
end

function Y = WaveletHandle(X, lvl, wname, xdim)
    [Y, ~] = wavedec2(reshape(X, xdim), lvl, wname);
    Y = Y(:);
end

function Y = WaveletHandleT(X, S, wname, xdim)
    Y = waverec2(X, S, wname);
    Y = reshape(Y, xdim);
end

function Y = Afun(X, flag, kernelF, S, K, xdim, indices, convolve_flag)
     if convolve_flag
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
     else  
        if strcmp(flag,'notransp')
            %Nontransposed 2DConvolution and Subsampling
            x = X(1:end-1);
            mu = X(end);
    %        Y = ifft2(kernelF .* fft2(reshape(x, xdim)));
            Y =  S*x;

    %         z = reshape(x, xdim);
    %         imagesc(z)
            %Y = Y(indices); + mu;
        else
            %Transposed 2DConvolution and Subsampling
            tmp = zeros(xdim);
            tmp(indices) = X;
    %         Y = (ifft2(conj(kernelF) .* fft2(tmp)));
            tmp = reshape(tmp, [], 1);
            Y =  K*X;
            %imagesc(Y)
            Y = [Y(:); sum(X(:))];
        end
     end
end

