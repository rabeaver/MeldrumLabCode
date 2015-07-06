function matsoft = soften(mat,dim)

for i = 2 : dim-1
    max = mat(i-1)*mat(i);
    for j = i : (i+1)
        if mat(j)*mat(j)>max
            max = mat(j)*mat(j);
        end
    end
    matsoft (i) = sqrt(max);
end
matsoft(1) = mat(1);
if(mat(2))^2 > (mat(1))^2
    matsoft(1) = mat(2);
end
matsoft(dim) = mat(dim);
if (mat(dim-1))^2 > (mat(dim)^2)
    matsoft(dim) = mat(dim-1);
end

            