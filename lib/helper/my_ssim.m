function result = my_ssim(X, Y)
if max(X(:)) > 10
    result = ssim(uint8(X), uint8(Y));
else
    result = ssim(X, Y);
end
end