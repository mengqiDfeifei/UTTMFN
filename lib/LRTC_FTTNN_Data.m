function [Lk, RSE_list] = LRTC_FTTNN_Data(D, mask, X, opt)
%% parameters
ranks=opt.ranks;
lambda=opt.lambda;
tol=opt.tol;
maxIter=opt.maxIter;
transform=opt.transform;
rho=opt.rho;
mu=opt.mu;
max_mu=opt.max_mu;
[n1, n2, n3] = size(D);
%% Initialization
Lk  = D;
Lk_transform = lineartransform(Lk, transform);
errList = zeros(maxIter, 1);
normT   = norm(D(:));
Y{1}=zeros(n1,ranks,n3);
for i = 1:n3
    [u,s,v] = svd(Lk_transform(:,:,i));
    Uk_transform(:,:,i) = u(:,1:ranks)*sqrt(s(1:ranks,1:ranks));
    Vk_transform(:,:,i) = sqrt(s(1:ranks,1:ranks))*v(:,1:ranks)';
end
Uk_hat_transform = Uk_transform;
Vk_hat_transform = Vk_transform;
Y{2}=zeros(ranks,n2,n3);
Y{3} = zeros(n1, n2, n3);
RSE_list = [];
%% Iteration
for Iter = 1: maxIter
    if mod(Iter, 50) == 0
        fprintf('DTTNN: iterations = %d  beta = %d  difference=%f time=%f\n', Iter, mu, errList(Iter-1),toc);
    end
    if isfield(opt, 'DEBUG')
        if opt.DEBUG
            time = toc;
            temp = [norm(Lk(:)-X(:), "fro")/ norm(X(:),"fro");time];
            
            RSE_list = [RSE_list, temp];
            fprintf('RSE: %f, time=%f\n', temp(1), temp(2));
        end
    end
    Lk_pre = Lk;
    for i = 1:n3
        % Update Uk and Vk
        Uk_transform(:,:,i) = (mu*Uk_hat_transform(:,:,i)+(mu*Lk_transform(:,:,i)+Y{3}(:,:,i))*Vk_transform(:,:,i)'-Y{1}(:,:,i))/((1+mu)*eye(ranks)+mu*Vk_transform(:,:,i)*Vk_transform(:,:,i)');
        Vk_transform(:,:,i) = inv((1+mu)*eye(ranks)+mu*Uk_transform(:,:,i)'*Uk_transform(:,:,i))*(mu*Vk_hat_transform(:,:,i)+Uk_transform(:,:,i)'*(mu*Lk_transform(:,:,i)+Y{3}(:,:,i))-Y{2}(:,:,i));
        % Update Uk_hat and Vk_hat
        Uk_hat_transform(:,:,i) = (mu*Uk_transform(:,:,i)+Y{1}(:,:,i))/(2/3+mu);
        Vk_hat_transform(:,:,i) = my_svt(Vk_transform(:,:,i)+Y{2}(:,:,i)/mu, 2/(3*mu));
        % Update Lk
        Lk_transform(:,:,i) = Uk_transform(:,:,i)*Vk_transform(:,:,i)-Y{3}(:,:,i)/mu;
        Y{1}(:,:,i) = Y{1}(:,:,i)+mu*(Uk_transform(:,:,i)-Uk_hat_transform(:,:,i));
        Y{2}(:,:,i) = Y{2}(:,:,i)+mu*(Vk_transform(:,:,i)-Vk_hat_transform(:,:,i));
    end
    temp = (lambda*D+mu*tprod(inverselineartransform(Uk_transform, transform), inverselineartransform(Vk_transform, transform), transform)-inverselineartransform(Y{3},transform))/(lambda+mu);
    Lk = inverselineartransform(Lk_transform, transform);
    Lk(mask) = temp(mask);
    Lk_transform = lineartransform(Lk, transform);
    for i=1:n3
        Y{3}(:,:,i) = Y{3}(:,:,i)+mu*(Lk_transform(:,:,i)-Uk_transform(:,:,i)*Vk_transform(:,:,i));
    end
    %% Stopping criterion
    errList(Iter) = norm(Lk(:)-Lk_pre(:))/normT;
    if errList(Iter) < tol
        break;
    else
        mu = min(mu * rho, max_mu);
    end
end
Lk(mask) = X(mask);
end