function y = rfft(x)
     fft_x = fft(x);
     y = fft_x(1:(floor(length(fft_x)/2)+1));
end