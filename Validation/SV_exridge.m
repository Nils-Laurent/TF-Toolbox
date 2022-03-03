function [] = SV_exridge()

%% MCS signal A
Lx = 1024;
Nfft = 1024;
% sigma_w = 0.05;
sigma_w = 0.03;

Tx = (0:Lx-1)/Lx;

a1 = exp(2*(1-Tx).^3 + 1.5*Tx.^4);
a2 = 1+ 5*Tx.^3 + 7*(1-Tx).^6;

phi1 = 50*Tx+30*Tx.^3-20*(1-Tx).^4;
phi2 = 340*Tx-2.*exp(-2*(Tx-0.2)).*sin(14*pi.*(Tx-0.2));

IF1 = 50+90*Tx.^2+80*(1-Tx).^3;
IF2 = 340+4*exp(-2*(Tx-0.2)).*sin(14*pi.*(Tx-0.2))-28*pi.*exp(-2*(Tx-0.2)).*cos(14*pi.*(Tx-0.2));

s1 = a1.*exp(2*pi*1i*(phi1));
s2 = a2.*exp(2*pi*1i*(phi2));

s = s1+s2;

[g, Lh] = gauss_win(Lx, sigma_w);

V = stft(s, Nfft, g);
[Cs, Es] = exridge(V, 10, 2, 10);

X_win = (Lh + 1):(Lx - Lh);

M1 = max(abs(Cs(1, X_win) - IF1(X_win)));
M2 = max(abs(Cs(2, X_win) - IF2(X_win)));

if M1 < 2
    fprintf("[OK] RD mode 1\n");
else
    error("error is too large : max|ridge - IF| = %d", M1);
end
if M2 < 18
    fprintf("[OK] RD mode 2\n");
else
    error("error is too large : max|ridge - IF| = %d", M2);
end

end