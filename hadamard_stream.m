n=hadamard(96);
nn=ones(96);
n=0.5*(n+nn);
n = padarray(n,[0 16],0,'both');

imagesc(n)