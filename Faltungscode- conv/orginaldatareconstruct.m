
addpath('E:\Faltungscode_ver310818\');
Masken = 'convolutionmasks bksize 5 full resolution';
load(Masken);
opt = struct;
opt.tol = 1e-18; %tolerance
opt.max_time = 10 * 60; %maximum calculation time in seconds
opt.max_iter = 5000; %maximimum iterations

% For flexbox
opt.method ='flexbox';%'flexbox';
opt.flexbox_internal_iter = 50;  %How many iterations before checkpoint -> callback
opt.flexbox_nnls = 1; %use nnls or nor
opt.flexbox_l1grad = 0;%; %1e5; %factor for ||grad(f)||_1 regularization (this is TV-1-Norm)
opt.flexbox_l2grad = 0; %factor for ||grad(f)||_2^2 regularization (TV-2-Norm)
opt.flexbox_wname = 'none'; %use some wavelet or 'none'
 convolve_flag = 1;
 A = zeros;
[x, mu] = reconstruct(kernel,A,transpose(A) ,indices(1:length(data)), data(1:length(data)), convolve_flag, opt);
   imagesc(x)