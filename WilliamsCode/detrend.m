function y = detrend(input_array)
    x = linspace(0,length(input_array)-1,length(input_array));
    m = (sum(x.*input_array)-sum(x)*sum(input_array)/length(x))/(sum(x.^2)-sum(x).^2/length(x));
    b = mean(input_array) - m*mean(x);
    y = input_array - (m*x + b);
end