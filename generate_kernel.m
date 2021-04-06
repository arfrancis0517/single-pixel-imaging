function matrix = genkernel_stream(probability, blockSize, dimH, dimW, method)
% Returns a dimH x dimW logical 0-1 matrix, where each block is equal to 1
% with the given probability.

    if strcmp(method, 'stream')
        littleDimH = ceil(dimH / blockSize);
        littleDimW = ceil(dimW / blockSize);
        littleMatrix = (rand(littleDimH, littleDimW) < probability);
        matrix = kron(littleMatrix, ones(blockSize, blockSize, 'double'));
        matrix = matrix(1:dimH, 1:dimW);
    elseif strcmp(method, 'circles')
        matrix = rand(dimH, dimW) < probability;
        circ = zeros(dimH, dimW);
        tmp = fspecial('disk', blockSize);
        circ(1:size(tmp,1), 1:size(tmp,2)) = fspecial('disk', blockSize) > 1e-12;
        circ = fftshift(circ);
        matrix = conv2(circ, matrix, 'same') > 0.5;
        
        
    elseif strcmp(method, 'triangle')
        
            littleDimH = ceil(dimH / blockSize);
        littleDimW = ceil(dimW / blockSize);
%         littleMatrix = (rand(littleDimH, littleDimW) < probability);
        L=rand(1,littleDimH*littleDimW)< probability;
        
        
        M=ones(blockSize, blockSize, 'double' );
        Matrix=cell(blockSize);
    
        n=randperm(littleDimH*littleDimW);
        
         j=1;
    
        Matrix_vert= L(1+littleDimH*(j-1))*rot90(tril(M), n(1+littleDimH*(j-1)));
        
        for i=2+littleDimH*(j-1):littleDimH+littleDimH*(j-1)
      
        t=L(i)*rot90(tril(M), n(i));
        Matrix_vert=vertcat(Matrix_vert,t);
        end 
        Matrix_horz = Matrix_vert;
        
        
        for j=2:littleDimW
    
        Matrix_vert= L(1+littleDimH*(j-1))*rot90(tril(M), n(1+littleDimH*(j-1)));
        
        for i=2+littleDimH*(j-1):littleDimH+littleDimH*(j-1)
     
        t=L(i)*rot90(tril(M), n(i));
        Matrix_vert=vertcat(Matrix_vert,t);
        end 
        Matrix_horz =horzcat(Matrix_horz,Matrix_vert);
        end
        matrix = Matrix_horz;
        
        
    elseif strcmp(method, 'sym_triangle')
        
        littleDimH = ceil(dimH / blockSize);
        littleDimW = ceil(dimW / blockSize);
        L=rand(1,littleDimH*littleDimW)< probability; 
        M=zeros(blockSize, blockSize, 'double' );
        Matrix=cell(blockSize);
        n=randperm(littleDimH*littleDimW);   
        r=randi([1 90], 1, littleDimH*littleDimW);
        
        
        
        [xx,yy] = meshgrid(1:blockSize,1:blockSize);
        theta = atan2d(xx,yy);  
        x=blockSize;
        
        M(atan2d(xx,yy) <=r(1) ) = 1;
        
        B = [fliplr(M) M];
        
        M = B(:,x/2+1:x+x/2);
        
        
        j=1;
            
        Matrix_vert= L(1)*rot90(M, n(1));
        
        for i=2+littleDimH*(j-1):littleDimH+littleDimH*(j-1)
        M=zeros(blockSize, blockSize, 'double' );
        M(atan2d(xx,yy) <= r(i) ) = 1;
     
        B = [fliplr(M) M];
        
        M = B(:,x/2+1:x+x/2);
      
        
        t=L(i)*rot90(M, n(i));
        Matrix_vert=vertcat(Matrix_vert,t);
              
        end 
        
        Matrix_horz = Matrix_vert;
        
        
         for j=2:littleDimW
        M=zeros(blockSize, blockSize, 'double' );
        M(atan2d(xx,yy) <= r(1+littleDimH*(j-1)) ) = 1;
        
        B = [fliplr(M) M];
        
        M = B(:,x/2+1:x+x/2);
    
        Matrix_vert= L(1+littleDimH*(j-1))*rot90(M, n(1+littleDimH*(j-1)));
        
        for i=2+littleDimH*(j-1):littleDimH+littleDimH*(j-1)
           % i=1;
        M=zeros(blockSize, blockSize, 'double' );
        M(atan2d(xx,yy) <= r(i) ) = 1;
        
        B = [fliplr(M) M];
        
        M = B(:,x/2+1:x+x/2);
        
        t=L(i)*rot90(M, n(i));
        Matrix_vert=vertcat(Matrix_vert,t);
        end 
        Matrix_horz =horzcat(Matrix_horz,Matrix_vert);
         end


          matrix = Matrix_horz;


        
        
    end
        
        
 
end




% function matrix = genkernel_hadamard(dim)
% % Returns a dimH x dimW logical 0-1 matrix, where each block is equal to 1
% % with the given probability.
% 
% matrix = hadamard (dim);
% 
% 
%     if strcmp(method, 'stream')
%         littleDimH = ceil(dimH / blockSize);
%         littleDimW = ceil(dimW / blockSize);
%         littleMatrix = (rand(littleDimH, littleDimW) < probability);
%         matrix = kron(littleMatrix, ones(blockSize, blockSize, 'double'));
%         matrix = matrix(1:dimH, 1:dimW);
%     elseif strcmp(method, 'circles')
%         matrix = rand(dimH, dimW) < probability;
%         circ = zeros(dimH, dimW);
%         tmp = fspecial('disk', blockSize);
%         circ(1:size(tmp,1), 1:size(tmp,2)) = fspecial('disk', blockSize) > 1e-12;
%         circ = fftshift(circ);
%         matrix = conv2(circ, matrix, 'same') > 0.5;
%     end
%         
%         
%  
% end