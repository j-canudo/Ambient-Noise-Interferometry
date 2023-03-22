function spc = whiten2(spc,fs,low,high)
    N = (length(spc)-1)*2;
    i1 = ceil(low/(fs/N));
    i2 = ceil(high/(fs/N));
    spc(1:i1) = cos(linspace(pi/2,pi,i1)).^2 .* exp(1i*angle(spc(1:i1)));
    spc(i1:i2) = exp(1i*angle(spc(i1:i2)));
    spc(i2:end) = cos(linspace(pi,pi/2,length(spc)-i2+1)).^2 .* exp(1i*angle(spc(i2:end)));
end