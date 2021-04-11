function masks = genmask_cconv2(kernel, targetDim, indices)
	%targetDim is the size of the input vector x
	%indices can be an integer vector with values in 1 ... prod(targetDim)
	
	%Pad or truncate the kernel
    if sum(targetDim ~= size(kernel)) > 0
        tmp = kernel;
        kernel = zeros(targetDim, 'like', tmp);
        kernel(1:min(targetDim(1),size(tmp,1)), 1:min(targetDim(2),size(tmp,2))) = tmp(1:min(targetDim(1),size(tmp,1)), 1:min(targetDim(2),size(tmp,2)));
    end
	
	%flip the kernel
	kernel = rot90(kernel,2);
    
	
	%Preallocate and take over the data type of kernel
	masks = zeros(targetDim(1), targetDim(2), length(indices), 'like', kernel);
	
	for k = 1:length(indices)
		ind = indices(k) - 1;
		i = mod(ind, targetDim(1));
		j = floor(ind / targetDim(1));
		
		%This is simply a shifted variant of the flipped kernel each time!
		masks(:,:,k) = circshift(kernel, [1+i, 1+j]);
        imagesc(masks)
	end
end

