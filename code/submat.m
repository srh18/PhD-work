function mat = submat(mat,i,val)
mat = mat(mat(:,i) ==val,:);
end
