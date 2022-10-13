function result = my_psnr(X, Y)
if max(X(:)) > 10
    result = psnr(uint8(X), uint8(Y));
else
    result = psnr(X, Y);
end
end