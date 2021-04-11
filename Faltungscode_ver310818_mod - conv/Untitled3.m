kernel = im2double(imread('Kernel16.tif'));
obj1= im2double(imread('USAF-3.tif'));
kernel_fft = fft2(kernel);
imagesc(obj1)
target_fft = fft2(obj1);

yr = ifft2(kernel_fft .* target_fft);
y = reshape(yr,1,[]);
h= lsqr(kernel,y(1:48)', 1e-20,100);
imagesc();