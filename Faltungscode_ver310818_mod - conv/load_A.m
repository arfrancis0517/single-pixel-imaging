
X = 64; %64
Y = 48; %48

M = 2048;

masks_dir = 'C:\Computational Imaging Simulations\masks\cs-comp 16\';
A=zeros(X*Y,M);

for i=1:M;  
    temp=imread([masks_dir num2str(i) '.png']);
    temp=reshape(temp,X*Y,1);
    
    A(:,i) = temp;
end
imshow (A);