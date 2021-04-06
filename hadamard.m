%save('C:\Users\admkRMC-Cosip\Desktop\Hadamard Sim\test.mat')
%{%}
clear all
close all
clc;tic
load ('C:\Users\admkRMC-Cosip\Desktop\Guo\09.01.2020\DLR_logo_90°_Hardamad_128_0.1_start_with_5000_pseudo_masks.mat','data');
N=128;

figure;
plot(data);

% k=1;
% for i=4:2:length(data)%6
%     y(k)=data(i);
%     k=k+1;
% end

k=1;
yplus=zeros(N*N,1);
yminus=zeros(N*N,1);
for i=2+5000:2:(length(data)-1)
    yplus(k)    = abs(data(i));
    yminus(k)   = abs(data(i+1));
    k=k+1;
end

y= -(yplus-yminus);


erg=ifwht(y);
res=reshape(erg,N,N);
figure;
imagesc(res), axis tight, colormap jet,axis equal;
figure;
plot(y);