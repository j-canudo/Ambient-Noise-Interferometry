function tr_out = fk_filt(tr_in,fs,dx,sgn,vel)
    if isempty(vel)
        cmin = 5;
        cmax = 50;
    else
        cmin = vel(1);
        cmax = vel(2);
    end
    if isempty(sgn)
        sgn = 'pos';
    end
    Nx = size(tr_in,1);
    Ns = size(tr_in,2);
    f0 = -fs/2:fs/Ns:fs/2-fs/Ns;
    k0 = -1/dx/2:1/dx/Nx:1/dx/2-1/dx/2/Nx;
    ft2 = fftshift(fft2(tr_in));
    [F,K] = meshgrid(f0,k0);
    C = F./K;
    filt = zeros(size(ft2));
    if isequal(sgn,'pos')
        filt((C > cmin) & (C < cmax)) = 1;
    else
        filt((C < -cmin) & (C > -cmax)) = 1;
    end
    filt = imgaussfilt(filt,3);
    ft2f = ft2.*filt;
    tr_out = real(ifft2(fftshift(ft2f)));
end