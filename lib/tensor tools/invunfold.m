function result = invunfold(X,n1,n2,n3)
result = zeros(n1,n2,n3);
for i = 1:n3
    temp = reshape(X(i,:),[n2,n1]);
    result(:,:,i) = temp';
end

end