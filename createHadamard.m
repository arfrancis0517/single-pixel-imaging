N=128;
k=1;
l=1;
for i=1:800
    % mask=Phi(i,:);
    delta=zeros(1,N*N);
    mask=zeros(1,N*N);
    delta(i)=1;
    aux=ifwht(delta);
    mask=1/2*(1+aux);
    mask=reshape(mask,N,N);
   
    MASK=imresize(mask,[1024 1024],'nearest');
    MASK=MASK(1:768,:);
    
    imwrite(MASK,['/Users/francis/Desktop/Hadamard/' num2str(k), '.png'])
    k=k+1;
    l=l+1
    imwrite(1-MASK,['/Users/francis/Desktop/Hadamard/' num2str(k), '.png'])
    k=k+1;
    l=l+1
end
% k=2;
% for i=174751:174751+174750-1
%     % mask=Phi(i,:);
%     delta=zeros(1,N*N);
%     mask=zeros(1,N*N);
%     delta(i)=1;
%     aux=ifwht(delta);
%     mask=1/2*(1+aux);
%     mask=reshape(mask,N,N);
%     MASK=zeros(768,1024);
%     MASK=imresize(mask,[1024 1024],'nearest');
%     MASK=MASK(1:768,:);
%     
%     imwrite(MASK,['C:\Users\admkRMC-Cosip\Desktop\Hadamard Sim\Streams\' int2str(N) '\2\' int2str(k),'.png'])
%     k=k+1;
%     l=l+1
% end
% k=2;
% for i=349501:349501+174750-1
%     % mask=Phi(i,:);
%     delta=zeros(1,N*N);
%     mask=zeros(1,N*N);
%     delta(i)=1;
%     aux=ifwht(delta);
%     mask=1/2*(1+aux);
%     mask=reshape(mask,N,N);
%     MASK=zeros(768,1024);
%     MASK=imresize(mask,[1024 1024],'nearest');
%     MASK=MASK(1:768,:);
%     
%     imwrite(MASK,['C:\Users\admkRMC-Cosip\Desktop\Hadamard Sim\Streams\' int2str(N) '\hilfs\' int2str(k),'.png'])
%     k=k+1;
%     l=l+1
% end
% 
