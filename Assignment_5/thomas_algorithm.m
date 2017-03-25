function [ output_args ] = thomas_algorithm( mat_A,b )
    
    l = length(mat_A);
    
    for i=1:l-1,
        fac = mat_A(i+1,i) / mat_A(i,i);
        mat_A(i+1,:) = mat_A(i+1,:) - fac*mat_A(i,:);
        b(i+1,:) = b(i+1,:) - fac*b(i,:);
    end
    for i=l:-1:2,
        fac = mat_A(i-1,i) / mat_A(i,i);
        mat_A(i-1,:) = mat_A(i-1,:) - fac*mat_A(i,:);
        b(i-1,:) = b(i-1,:) - fac*b(i,:);
    end
    
    output_args = b ./ diag(mat_A);

end