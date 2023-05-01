%% Tensor Contraction for A_ijk * B_k

%% Function
function result = TenCont(A,B,dim)
    % Matrix Size
    sizeA = size(A);
    sizeB = size(B);
    % Error of dimension
    if size(dim) > 2
        disp('Wrong Assignment of Dimension in TenCont');
    end
    % Size of the Target Column
    sizetar = max(sizeB);
    % Alert of Illegal Operation
    if sum(sizeB) > sizetar + 1
        dispstat(sprintf('Illegal Tensor Contraction'),'keepthis','timestamp');
    end
    % Erasing the Contribution of the Target Column
    sizeA(dim(1)) = [];
    % Reshaping Matrix and Contraction
    result = reshape(reshape(A,[prod(sizeA),sizetar])*B,sizeA);
end