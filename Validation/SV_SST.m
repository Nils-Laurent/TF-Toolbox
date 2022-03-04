function [] = SV_SST()

%% s1
Lx = 1024;
Nfft = 1024;
Tx = (0:(Lx-1))/Lx;
s = exp(2i*pi*(200*Tx));
sigma_w = 0.04;

[~, SSTN, ~] = sstn(s, sigma_w, Nfft);

[~, SST2, ~] = sst2(s, sigma_w, Nfft);

if max(max(abs(SSTN.d1 - SST2.d1))) > 0 ||...
    max(max(abs(SSTN.d2 - SST2.d2))) > 0
    error("SSTs are not equal");
else
    fprintf("[OK] sst2.d1 = sstn.d1 (same for d2)\n");
end

%% s2
Lx = 2048;
Nfft = 512;
sigma_w = 0.03;

Tx = (0:Lx-1)/Lx;

% a1 = exp(2*(1-Tx).^3 + 1.5*Tx.^4);
% a2 = 1+ 5*Tx.^3 + 7*(1-Tx).^6;

phi1 = 50*Tx+30*Tx.^3-20*(1-Tx).^4;
phi2 = 340*Tx-2.*exp(-2*(Tx-0.2)).*sin(14*pi.*(Tx-0.2));

sA = exp(2*pi*1i*(phi1));
sB = exp(2*pi*1i*(phi2));

s2 = sA+sB;

[~, TN, ~] = sstn(s2, sigma_w, Nfft);
[~, Lh] = gauss_win(Lx, sigma_w);
X_win = (Lh+1):(Lx-Lh);

r3 = zeros(1, Lx);
r4 = r3;
for n=1:Lx
    % g(0) = 1, no need to divide by g(0)
    [~, ii3] = sort(abs(TN.d3(:, n)), 'desc');
    r3(n) = sum(TN.d3(ii3(1:2), n))/Nfft;
    [~, ii4] = sort(abs(TN.d4(:, n)), 'desc');
    r4(n) = sum(TN.d4(ii4(1:2), n))/Nfft;
end

EM3 = max(abs(r3(X_win) - s2(X_win)));
EM4 = max(abs(r4(X_win) - s2(X_win)));

if (EM3 > 1) || (EM4 > 1)
    error("SST reconstruction error too big");
else
    fprintf("[OK] sstn\n");
end

end

