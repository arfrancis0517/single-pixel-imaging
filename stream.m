function matrix=stream(ind,p,Psize)

matrix=zeros(1024,768);
for i=1:Psize:1024
    for j=1:Psize:768
                aux = binornd(1,p);
        for k=i:i+Psize-1
            for l=j:j+Psize-1
                matrix(k,l)=aux;
            end
        end
        
    end
end
h=figure;
%matrix=matrix./max(max(matrix));
matrix=matrix';
matrix=logical(matrix);
set(gca,'position',[0 0 1 1],'units','normalized')
set(gcf,'units','pixel','PaperPosition',[0 0 1024 768]);
%imshow(matrix)[0 0 1024 768]
grid off
axis off
colormap gray

name=['C:\Users\admkRMC-Cosip\Desktop\Guo\Convolutional_Blksize_16_100k_sym_triangle_p_0.5\' num2str(ind)];
%print(h,'-dpng',name,'-r1')
close all
imwrite(matrix,[name,'.png'])

%saveSameSize(gcf, 'format', 'bmp'); 