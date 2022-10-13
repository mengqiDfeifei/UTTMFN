function Y=unfold(X,n)
[n1,n2,n3]=size(X);
if n==1
    
elseif n==2
    
elseif n==3
    Y=zeros(n3,n1*n2);
    k=1;
    for i=1:n1
        for j=1:n2
            Y(:,k)=X(i,j,:);
            k=k+1;
        end
    end
end
end