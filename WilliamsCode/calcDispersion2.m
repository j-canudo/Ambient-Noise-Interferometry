function [frq,vel,disp] = calcDispersion2(tr_,fs,dx,v_min,v_max,v_step,f_min,f_max)
    vel = v_min:v_step:v_max;
    Nv = length(vel);
    Nx = size(tr_,1);
    Ns = size(tr_,2);
    for i=1:Nx
        sp_(i,:) = rfft(tr_(i,:));
    end
    frq = linspace(0,fs/2,floor(Ns/2)+1);
    sp_shift = zeros(Nx,floor(Ns/2)+1,Nv);
    disp = zeros(Nv,floor(Ns/2)+1);
    for iv=1:Nv
        for ix=1:size(tr_,1)
            sp_shift(ix,:,iv) = sp_(ix,:) .* exp(2*1i*pi*frq*(ix-1)*dx/vel(iv));
        end
        disp(iv,:) = mean(real(sp_shift(:,:,iv)),1);
    end
end