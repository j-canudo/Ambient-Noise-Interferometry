function y = irfft(x)
    if mod(length(x),2) == 0
        even = true;
    else
        even = false;
    end
    if (even)
        n = 2 * (length(x) - 1 );
        s = length(x) - 1;
    else
        n = 2 * (length(x) - 1 );
        s = length(x);
    end
    xn = zeros(1,n);
    xn(1:length(x)) = x;
    xn(length(x):n) = conj(x(s:-1:2));
    y  = real(ifft(xn));
end