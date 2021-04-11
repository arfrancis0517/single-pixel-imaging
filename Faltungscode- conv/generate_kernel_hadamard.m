
function matrix = genkernel_hadamard(dim,scale, height, width)
% Returns a dimH x dimW logical 0-1 matrix, where each block is equal to 1
% with the given probability.

matrix = hadamard (dim,scale);
%matrix([769:1024],:)=[];
%matrix([768:1024],:)=[];
%matrix = matrix(1:1024,769:1024)

%matrix = matrix + 1; making Hadamard zeros and one rather than -1 and 1s
%matrix = logical(matrix);
%imagesc(matrix)


        
        
 
end