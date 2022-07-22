function [m1,m2,m3,mat1,mat2,mat3] = mergemats(del)
mat = load_good_mat(2,1);
mat = submat(mat,2,del);
mat1 = submat(load_good_mat(3,1),2,del);
mat2 = submat(load_good_mat(4,1),2,del);
mat3 = submat(load_good_mat(5,1),2,del);
m1 = submat(mat,3,1); 
m2 = submat(mat,3,2); 
m3 = submat(mat,3,3); 
end


function mat = submat(mat,i,val)
mat = mat(mat(:,i) ==val,:);
end
