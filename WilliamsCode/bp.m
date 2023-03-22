function y = bp(x,freqmin,freqmax,fs)
    [b,a] = butter(4,[freqmin freqmax]/(fs/2),'bandpass');
    y = filtfilt(b,a,x);
end