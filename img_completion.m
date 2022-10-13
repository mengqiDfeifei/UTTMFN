clear;clc;close all;
addpath(genpath(cd))
%%
sr_all=[0.05,0.15,0.3];
r_all=[30,50,100];
X = double(imread("re1.jpg"));
[n1,n2,n3]=size(X);
for sr=sr_all
    for r = r_all(sr==sr_all)
        %% Input data
        disp([num2str((1-sr)*100) '% pixels are missing']);
        temp=double(rand(n1, n2, n3) < sr);
        omega = logical(temp);
        M = zeros(n1, n2, n3);
        M(omega) = X(omega);
        %% parameter settings
        opts.mu=1e-4;
        opts.tol=1e-4;
        opts.rho = 1.05;
        opts.DEBUG = 0;
        opts.max_mu = 1e20;
        opts.lambda = 1/sqrt((max(n1, n2)*n3));
        opts.maxIter = 500;
        M = interpolation(M, omega);
        [transform.L,~,~]=svds(unfold(M,3),n3);
        transform.L=transform.L';
        transform.l=1;
        transform.inverseL=inv(transform.L);
%         transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
%         transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
        opts.transform=transform;
        opts.ranks=r;
        %% UTTMFN_1_2
        tic
        UTTMFN_1_2 = LRTC_DTTNN_Data(M, omega, X, opts);
        time = toc;
        UTTMFN_1_2_PSNR = my_psnr(UTTMFN_1_2, X)
        UTTMFN_1_2_SSIM = my_ssim(UTTMFN_1_2, X)
        disp(['UTTMFN-1/2',' time:', num2str(time),' psnr:', num2str(UTTMFN_1_2_PSNR), ' ssim:', num2str(UTTMFN_1_2_SSIM)])
        %% UTTMFN_2_3
        tic
        UTTMFN_2_3 = LRTC_FTTNN_Data(M, omega, X, opts);
        time = toc;
        UTTMFN_2_3_PSNR = my_psnr(UTTMFN_2_3, X)
        UTTMFN_2_3_SSIM = my_ssim(UTTMFN_2_3, X)
        disp(['UTTMFN-1/2',' time:', num2str(time),' psnr:', num2str(UTTMFN_2_3_PSNR), ' ssim:', num2str(UTTMFN_2_3_SSIM)])
        %% UTTMFN_1_3
        tic
        UTTMFN_1_3 = LRTC_Tr_TTNN_Data(M, omega, X, opts);
        time = toc;
        UTTMFN_1_3_PSNR = my_psnr(UTTMFN_1_3, X)
        UTTMFN_1_3_SSIM = my_ssim(UTTMFN_1_3, X)
        disp(['UTTMFN-1/2',' time:', num2str(time),' psnr:', num2str(UTTMFN_1_3_PSNR), ' ssim:', num2str(UTTMFN_1_3_SSIM)])
    end
end