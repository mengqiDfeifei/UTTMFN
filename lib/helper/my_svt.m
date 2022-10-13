function result=my_svt(temp,lambda)

[temp_U, temp_S, temp_V] = svd(temp, 'econ');
temp_S = diag(temp_S);
svp = length(find(temp_S > lambda)); 
if svp>=1
    temp_S = temp_S(1:svp) - lambda;
else
    svp   = 1;
    temp_S = 0;
end
result = temp_U(:,1:svp)*diag(temp_S(1:svp))*temp_V(:,1:svp)';
end