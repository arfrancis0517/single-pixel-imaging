
% clear all
% close all
% clc;tic
%load ('C:\Users\admkRMC-Cosip\Desktop\Guo\08.01.2020\DLR_logo_90°_Pseudo_16_mesurements_7500_0.05_white_ref.mat','data');
load ('C:\Users\admkRMC-Cosip\Desktop\Guo\matfile_maks\mask_W_Str_P_0.5_Bl_16.mat','m');

warning off
Mx=64;%32;%64
My=48;%24;%48
%M=20000;

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

figure;
plot(data,'.')

figure;
plot(y,'.')

a = median(white(:))./median(y(:));
Sig = mad(y(:)*a);
ref = mad(white);
SNR = Sig/ref


auxi_w = white;

% drift correction 
for i=1:length(white)
y(i)=y(i).*auxi_w(5000)./auxi_w(i);
white(i)=auxi_w(i).*auxi_w(5000)./auxi_w(i);
% black(i)=black(i).*auxi_w(2000)./auxi_w(i);
end

figure;
plot(y,'.')

figure
plot(white,'.','Color',[0.831 0.816 0.784])
hold on
plot(y,'.','Color','blue')
%plot(black,'.','Color','black')
hold off



a = median(white(2000:length(y)))./median(y(2000:length(y)));
Sig = mad(y(2000:length(y))*a);
ref = mad(white(2000:length(y)));
SNR = Sig/ref

j=1;
k=1;
for i=j:length(y)
 mask=imresize(m{i},1,'nearest');
 Phi(k,:)=  - reshape(mask,1,Mx*My);
 Phi=double(Phi);
 Phi(k,:) = white(k) + Phi(k,:);
 %Phi(i,:) = white(i) + white(i)*((black(i)/white(i)).*(Phi(i,:)));
 k=k+1;
end
% figure
% imagesc(Masken{4})

a = median(white(2500:length(y)))./median(y(2500:length(y)));
Sig = mad(y(2500:length(y))*a);
ref = mad(white(2500:length(y)));
SNR = Sig/ref


x=nnls(Phi(2500:length(y),:),y(2500:length(y))'); 

res=reshape(x,My,Mx);
% h=fspecial('gaussian', 2);
% figure;
% imagesc(h);
% n=imfilter(res,h);
% imagesc(n)
figure
imagesc(res./max(max(res))),axis tight, colormap jet,axis equal
toc;