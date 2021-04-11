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