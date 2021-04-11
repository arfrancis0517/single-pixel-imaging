clc;
clear all; 
close all;

% figure
% kernel=generate_kernel(0.0001,32,768,1024,'circles');
n=10000;

kernel=generate_kernel(0.5,32,768,1024,'stream');

% figure
% imagesc(kernel)
% drawnow

indices = randsample(768*1024, n);

for i=1:n
%figure
matrix = generate_mask(kernel,[768,1024],indices(i));
set(gca,'position',[0 0 1 1],'units','normalized')
set(gcf,'units','pixel','PaperPosition',[0 0 1024 768]);
%imshow(matrix)[0 0 1024 768]
grid off
axis off
colormap gray

name=['C:\Users\admkRMC-Cosip\Desktop\Guo\Convolutional_Blksize_32_100k\' num2str(i)];
%print(h,'-dpng',name,'-r1')
close all
imwrite(matrix,[name,'.png'])
end

