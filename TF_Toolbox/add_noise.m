function [y, noise, Scale] = add_noise(x, SNR_in, varargin)
%ADD_NOISE add noise to x such that the input SNR is SNR_in
%   [y, noise] = add_noise(x, SNR_in)
%   [y, noise] = add_noise(x, SNR_in, n)
%
% INPUTS:
%   x        : signal
%   SNR_in   : input SNR
%
%   n        : noise
%              default is a real randn if x is real
%                         or if parameter 'real' is true.
%                      is a complex randn if x has non-zero
%                         imaginary part or if parameter 'real' is true.
%
%   --  List of possible name-value pair argument
%   'real'   : no effect if n is given,
%              noise is real gaussian noise (see randn).
%   'imag'   : no effect if n is given,
%              noise is gaussian and has a nonzero imaginary part.
%
% OUTPUTS:
%   --  Output
%   y        : noisy signal.
%   noise    : noise added to the signal.

L = length(x);
defaultN = [];

p = inputParser;
addRequired(p, 'x');
addRequired(p, 'SNR_in');
addOptional(p, 'n', defaultN);
addParameter(p, 'real', false);
addParameter(p, 'imag', false);

parse(p, x, SNR_in, varargin{:});
r = p.Results;

if (r.imag == true && r.real == true)
  error("the noise connot be real and have nonzero imaginary part.");
end

if (sum(size(r.n)) > 0)
  %% noise is already provided
  % r.n = r.n
elseif (r.imag == true)
  %% imag flag is set
  r.n = randn(1, L) + 1i*randn(1, L);
elseif (r.real == true) || (sum(abs(imag(x))) == 0)
  %% real flag is set or signal real
  r.n = randn(1, L);
else
  r.n = randn(1, L) + 1i*randn(1, L);
end


x = x(:);
r.n = r.n(:);

if SNR_in == inf
    y = x(:);
    noise = zeros(size(y));
    return;
end

if SNR_in == -inf
    y = r.n(:);
    noise = y;
    return;
end

Scale = sqrt(sum(abs(x).^2)/sum(abs(r.n).^2)*10^(-SNR_in/10));
y = x + Scale*r.n;
noise = r.n;
end

