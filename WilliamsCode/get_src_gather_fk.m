function trxc = get_src_gather_fk(data, src, rec_arr, nnx, nns, nnw, nwn, nov, fmin, fmax, fs, dx, sgn,vel)
    [~,src_idx] = min(abs(rec_arr-src));
    nc = length(rec_arr);
    trxc = zeros(nc,nns);
    spxc = zeros(nc,nnw);
    for n=1:nwn
        tr = data(rec_arr(1):rec_arr(length(rec_arr)),(n-1)*nov+1:((n-1)*nov + nns));
        for ic=1:nc
            tr(ic,:) = detrend(tr(ic,:));
            tr(ic,:) = bp(tr(ic,:),fmin,fmax,fs);
        end
        tr = fk_filt(tr,fs,dx,sgn,vel);
        sp = zeros(nc,nnw);
        for ic=1:nc
            sp(ic,:) = rfft(tr(ic,:));
            sp(ic,:) = whiten2(sp(ic,:),fs,fmin,fmax);
        end
        for ic=1:nc
            spxc(ic,:) = spxc(ic,:) + conj(sp(ic,:)) .* sp(src_idx,:);
        end
    end
    for ic=1:nc
        trxc(ic,:) = irfft(spxc(ic,:));
    end
end