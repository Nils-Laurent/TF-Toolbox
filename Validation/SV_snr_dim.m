function [] = SV_snr_dim()

N = [8, 5, 512, 8];
dim = 3;

y = 30*randn(1, N(dim));
x = (1:512) + y;

M = zeros(N);

for n1=1:N(1)
    for n2=1:N(2)
        for n4=1:N(4)
            M(n1, n2, :, n4) = x + 5*randn(1, N(dim));
        end
    end
end

Y = snr_dim(x, M, dim);

for n1=1:N(1)
    for n2=1:N(2)
        for n4=1:N(4)
            m_ = squeeze(M(n1, n2, :, n4));
            m_ = m_(:).';
            if abs(Y(n1, n2, n4) - snr(x, m_ - x)) > 0
                error("snr_dim : snrs are not equal");
            end
        end
    end
end

fprintf("[OK] snr_dim\n");

end