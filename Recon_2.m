% Variables provided by Labview
clc;
clear all; 
close all;

load ('C:\Users\admkRMC-Cosip\Desktop\Guo\05.02.2020\1.mat','data');
load ('C:\Users\admkRMC-Cosip\Desktop\Guo\matfile_maks\Convolutional_Blksize_32_50k_triangle_p_0.5.mat');

Mx=32
My=24

k=1;
for i=3:2:length(data)%6
    white(k)=data(i);
    k=k+1;
end
k=1;
for i=4:2:length(data)%6
    y(k)=data(i);
    k=k+1;
end

% k=1;
% for i=3:2:length(data)%6
%     ind(k)=indices(i);
%     k=k+1;
% end



figure;
plot(data,'.')

auxi_w = white;

% drift correction 
for i=1:length(white)
y(i)=y(i).*auxi_w(5000)./auxi_w(i);
white(i)=auxi_w(i).*auxi_w(5000)./auxi_w(i);
% black(i)=black(i).*auxi_w(2000)./auxi_w(i);
end

% matrix = generate_mask(kernel,[768,1024],ind(1));
% 
% figure
% imagesc(matrix)


figure
plot(white,'.','Color',[0.831 0.816 0.784])
hold on
plot(y,'.','Color','blue')
%plot(black,'.','Color','black')
hold off


    median(y)
	opt = struct;
	opt.tol = 1e-15; %tolerance
    opt.max_time = 10 * 600; %maximum calculation time in seconds
	opt.max_iter = 1000000; %maximimum iterations
    %checkLength = length(data)-startMasks;%-startMasks

	opt.method = 'lsqr';
	opt.flexbox_internal_iter = 100;  %How many iterations before checkpoint -> callback
	opt.flexbox_nnls = 0; %use nnls or not
	opt.flexbox_l1grad = 1e3; %factor for ||grad(f)||_1 regularization (this is TV-1-Norm)
    opt.flexbox_l2grad = 0.00; %factor for ||grad(f)||_2^2 regularization (TV-2-Norm)
	opt.flexbox_wname = 'none'; %use some wavelet or 'none'
% 	kernel = imresize(kernel,4,'nearest');
	[x, mu] = reconstruct(kernel, indices(2000:length(y)), y(2000:length(y))', opt); %[x, mu] = reconstruct(kernel, indices(1:length(data)), ym', opt);
	%imagesc(z);
%     [x, mu] = reconstruct(kernel, indices(2000:19300), y(2000:19300)', opt); %[x, mu] = reconstruct(kernel, indices(1:length(data)), ym', opt);
	res=x; %
%     z=reshape(x,768,1024);
figure
imagesc(-res./max(max(res))),axis tight, colormap jet,axis equal
% toc;

