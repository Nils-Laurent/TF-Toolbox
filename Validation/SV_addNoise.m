function [] = SV_addNoise()

SNRs = -10:1:20;
M = 30;

th = 1E-12;

Lx = 512;
x = ones(Lx, 1);

for snr_in = SNRs
    for m=1:M
        n = randn(1, Lx);
        y = add_noise(x, n, snr_in);
        snr_out = snr(x, y - x);
        if abs(snr_in - snr_out) > th
            error("snr diff too big\n");
        end
    end
end

fprintf("[OK] add noise\n");

end

