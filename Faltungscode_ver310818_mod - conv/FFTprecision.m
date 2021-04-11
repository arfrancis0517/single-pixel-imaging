%fftw('dwisdom',[]);
%method = fftw ('dwisdom')
kernel=(imread(['Kernel16.tif']));
kernel = (double(kernel));

L= length(kernel(:));
%NFFT = 2^nextpower2(L);
A= fft2(kernel).*fft2(x);
F = ifft2((A));
imagesc(F)
%imshow(F)

e =(abs(kernel(2,:)- F(2,:)));
plot(e)
class(kernel)
isequal(kernel, F)
subplot(2,1,1);
plot(kernel(2,:));
subplot(2,1,2);
plot(F(2,:));  

p=conv2(x,kernel,'same')