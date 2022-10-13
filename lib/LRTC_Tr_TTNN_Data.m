function [Lk, RSE_list] = LRTC_Tr_TTNN_Data(D, mask, X, opt)
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
    U1k_transform(:,:,i) = u(:,1:ranks);
    U2k_transform(:,:,i) = s(1:ranks,1:ranks);
    U3k_transform(:,:,i) = v(:,1:ranks)';
end
U1k_hat_transform = U1k_transform;
U2k_hat_transform = U2k_transform;
U3k_hat_transform = U3k_transform;
Y{2} = zeros(ranks, ranks, n3);
Y{3} = zeros(ranks,n2,n3);
Y{4} = zeros(n1, n2, n3);
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
        U1k_transform(:,:,i) = (mu*U1k_hat_transform(:,:,i)+(mu*Lk_transform(:,:,i)+Y{4}(:,:,i))*(U2k_transform(:,:,i)*U3k_transform(:,:,i))'-Y{1}(:,:,i))/((1+mu)*eye(ranks)+mu*(U2k_transform(:,:,i)*U3k_transform(:,:,i))*(U2k_transform(:,:,i)*U3k_transform(:,:,i))');
        A = U1k_transform(:,:,i)'*U1k_transform(:,:,i);
        B = U3k_transform(:,:,i)*U3k_transform(:,:,i)';
        C = U1k_transform(:,:,i)'*(Lk_transform(:,:,i)+Y{4}(:,:,i)/mu)*U3k_transform(:,:,i)'+U2k_hat_transform(:,:,i)-Y{2}(:,:,i)/mu;
        vecx = (kron(B',A)+kron(eye(ranks)',eye(ranks)))\C(:);
        U2k_transform(:,:,i) = reshape(vecx, ranks, ranks);
        U3k_transform(:,:,i) = pinv((1+mu)*eye(ranks)+mu*(U1k_transform(:,:,i)*U2k_transform(:,:,i))'*(U1k_transform(:,:,i)*U2k_transform(:,:,i)))*(mu*U3k_hat_transform(:,:,i)+(U1k_transform(:,:,i)*U2k_transform(:,:,i))'*(mu*Lk_transform(:,:,i)+Y{4}(:,:,i))-Y{3}(:,:,i));
        % Update Uk_hat and Vk_hat
        U1k_hat_transform(:,:,i) = my_svt(U1k_transform(:,:,i)+Y{1}(:,:,i)/mu, 1/(3*mu));
        U2k_hat_transform(:,:,i) = my_svt(U2k_transform(:,:,i)+Y{2}(:,:,i)/mu, 1/(3*mu));
        U3k_hat_transform(:,:,i) = my_svt(U3k_transform(:,:,i)+Y{3}(:,:,i)/mu, 1/(3*mu));
        % Update Lk
        Lk_transform(:,:,i) = U1k_transform(:,:,i)*U2k_transform(:,:,i)*U3k_transform(:,:,i)-Y{4}(:,:,i)/mu;
        Y{1}(:,:,i) = Y{1}(:,:,i)+mu*(U1k_transform(:,:,i)-U1k_hat_transform(:,:,i));
        Y{2}(:,:,i) = Y{2}(:,:,i)+mu*(U2k_transform(:,:,i)-U2k_hat_transform(:,:,i));
        Y{3}(:,:,i) = Y{3}(:,:,i)+mu*(U3k_transform(:,:,i)-U3k_hat_transform(:,:,i));
    end
    temp = tprod(inverselineartransform(U1k_transform, transform),inverselineartransform(U2k_transform, transform),transform);
    temp = tprod(temp, inverselineartransform(U3k_transform, transform), transform);
    temp = (lambda*D+mu*temp-inverselineartransform(Y{4},transform))/(lambda+mu);
    Lk = inverselineartransform(Lk_transform, transform);
    Lk(mask) = temp(mask);
    Lk_transform = lineartransform(Lk, transform);
    for i=1:n3
        Y{4}(:,:,i) = Y{4}(:,:,i)+mu*(Lk_transform(:,:,i)-U1k_transform(:,:,i)*U2k_transform(:,:,i)*U3k_transform(:,:,i));
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