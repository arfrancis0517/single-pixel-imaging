% Variables provided by Labview
clc;
clear all; 
close all;

% ON=1;
% state = 0;
% name = 'C:\Users\SvrenA\Desktop\convolutionmasks bksize 5 full resolution.mat';% data
% init = 1;
% Masken = 'convolutionmasks bksize 5 full resolution';
% ConvInd = 1;
% startRef = 0.0;
% startMasks = 0.0;% changed to zero, the dimensions were not matching
% Mx = 1024;
% My = 768;
% data = data * 10;
% 
% addpath('C:\Users\SvrenA\Desktop\convolutionmasks bksize 5 full resolution')
% load(Masken)
% save(name)
% M=length(data)/2-100;
% Commented only for convolutional
% k=1;
% for i=startRef:2:length(data)
%     ref(k)=data(i);
%     k=k+1;
% end
% k=1;
% for i=startMasks:2:length(data)
%     y(k)=data(i);
%     k=k+1;
% end
% 
%  
% Driftkorrektur !!! invertiert weiss = hoch, schwarz=niedrig0
% if state == 1
% 	for i=1:M
% 		y(i)=y(i).*ref(M)./ref(i);
% 		ref(i)=ref(i).*ref(M)./ref(i);
% 	end
% end
% 
% if ConvInd==0
% 
% 	schwarz(1)=data(1);
% 	for i=1:M
% 		Phi(i,:)=reshape(Masken{i},1,Mx*My);
% 		Phi=double(Phi);
% 		Phi(i,:) =schwarz(1) +schwarz(1)*(ref(1)/schwarz(1)-1)*Phi(i,:);
% 	end
% 	x=zeros(My*Mx,1);
% 	x=nnls(Phi(init:M-100,:),y(init:M-100)'); 
% 	res=reshape(x,My,Mx);
% 	
% else
% 	ym=zeros(1);
% 	s4=3;
% 	 
% 	
% 	for i=1:(length(data)- startMasks)
% 		if (state==4)
% 			startMasks=startMasks+2;
% 			ym(i)=data(startMasks);
% 			state=2;        
% 		else
% 			startMasks=startMasks+1;
% 			ym(i)=data(startMasks);        
% 			%state=state+1;        
% 		end
%     
%     end
%     
%     	
% 	for i=1:(length(data)- startMasks)
%         ym(i)=data(i+startMasks);
%     end
%   for i = 1:10000:400000
%    % data1 = data(200000:400000);
  

load ('C:\Users\admkRMC-Cosip\Desktop\Guo\28.01.2020\1.mat','data');
load ('C:\Users\admkRMC-Cosip\Desktop\Guo\matfile_maks\Convolutional_Blksize_32_100k.mat');

Mx=32;%64
My=24;%48

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

k=1;
for i=4:2:length(data)%6
    ind(k)=indices(i);
    k=k+1;
end



figure;
plot(data,'.')

auxi_w = white;

% drift correction 
for i=1:length(white)
y(i)=y(i).*auxi_w(5000)./auxi_w(i);
white(i)=auxi_w(i).*auxi_w(5000)./auxi_w(i);
% black(i)=black(i).*auxi_w(2000)./auxi_w(i);
end




figure
plot(white,'.','Color',[0.831 0.816 0.784])
hold on
plot(y,'.','Color','blue')
%plot(black,'.','Color','black')
hold off

   
	%Result = 'C:\Users\SvrenA\Desktop\Computational imaging\Masks verify rev 17 39 26_03_19\';
    median(y)
	opt = struct;
	opt.tol = 1e-15; %tolerance
    opt.max_time = 10 * 600; %maximum calculation time in seconds
	opt.max_iter = 1000000; %maximimum iterations
    %checkLength = length(data)-startMasks;%-startMasks

	opt.method = 'lsqr';
	opt.flexbox_internal_iter = 100;  %How many iterations before checkpoint -> callback
	opt.flexbox_nnls = 1; %use nnls or not
	opt.flexbox_l1grad = 1e3; %factor for ||grad(f)||_1 regularization (this is TV-1-Norm)
    opt.flexbox_l2grad = 0.00; %factor for ||grad(f)||_2^2 regularization (TV-2-Norm)
	opt.flexbox_wname = 'none'; %use some wavelet or 'none'
% 	kernel = imresize(kernel,4,'nearest');
	[x, mu] = reconstruct(kernel, ind, y', opt); %[x, mu] = reconstruct(kernel, indices(1:length(data)), ym', opt);
	imagesc(x);
    %saveas(s, fullfile(Result, [int2str(i),'.png']))
	res=reshape(x,My,Mx);
%end
% length(indices)
% 
% save(name)
