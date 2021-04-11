


kernel=imread('Kernel4x6.tif');
kernel=im2double(kernel);
imagesc(kernel)
height = 4;%for Bs 16   768;
width = 6;% for Bs 16  1024;
nmasks = 20;

%Get Target ready
target = mat2gray(imread('tar4x6.tif'));
imagesc(target)
%target = logical(target);
%target = phantom(1024);
%target = imresize(target, [height, width], 'Method', 'nearest');
%target = ones (height, width);


indices = randsample(height*width, nmasks);

kernel_fft = fft2(kernel);

target_fft = fft2(target);

yr = ifft2(kernel_fft .* target_fft);

imagesc(y)

y = yr(indices);

 imagesc(y)

subplot(3,1,1)
imagesc(kernel)
title('Kernel')
axis image

subplot(3,1,2)
imagesc(target)
title('Target')
axis image

subplot(3,1,3)
imagesc(yr)
title('Y')
axis image
