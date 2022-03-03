function [] = SV_downsamp()
%SV_DOWNSAMP short valid of the downsampled STFT
sigma_w = 0.04;

L = 4096;
t = (0:L-1)/L;
a  = 2;
s1 = a.*exp(2*pi*1i*(1000*t+60*cos(3*pi*t)));
s2 = a.*exp(2*pi*1i*(400*t+30*cos(3*pi*t)));
s  = s1+s2;
s = s(:);
Nfft = 512;

g = gauss_win(L, sigma_w);

cas = 1;
for downsamp=[1 4 8 16 32 64]
    s_vec = [0 1 2 10 20 40];
    for shift=s_vec
        if shift >= downsamp
            continue;
        end

        [tfr,~] = stft(s,Nfft,g,'cas',cas,'down',downsamp,'shift',shift);
        [x] = istft(tfr, g,'cas',cas,'len',length(s),'shift',shift);

        snr_out = snr(s, x - s);
        if (snr_out > 250)
            fprintf("[OK] snr_out=%f\n", snr_out);
        else
            error("snr_out=%f\n", snr_out);
        end
    end
end

downsamp = 32;
shift = 20;
cas = 2;
[tfr,~] = stft(s,Nfft,g,'cas',cas,'down',downsamp,'shift',shift);
[x] = istft(tfr, g,'cas',cas,'len',length(s),'shift',shift);

snr_out = snr(s, x - s);
if (snr_out > 250)
    fprintf("[OK] snr_out=%f\n", snr_out);
else
    error("snr_out=%f\n", snr_out);
end

cas = 3;
[tfr,~] = stft(s,Nfft,g,'cas',cas,'down',downsamp,'shift',shift);
[x] = istft(tfr, g,'cas',cas,'len',length(s),'shift',shift);

snr_out = snr(s, x - s);
if (snr_out > 250)
    fprintf("[OK] snr_out=%f\n", snr_out);
else
    error("snr_out=%f\n", snr_out);
end

end

