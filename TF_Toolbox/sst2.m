function [STFT,SST,omega] = sst2(s,sigma,Nfft,gamma)
%SST2 computes the STFT of a signal and different versions of synchrosqueezing
%   [STFT,SST,omega] = SST2(s,sigma,Nfft)
%   [STFT,SST,omega] = SST2(s,sigma,Nfft,gamma)
%
% INPUTS:
%   s     : Real or complex signal.
%   sigma : The parameter sigma in the definition of the Gaussian window.
%   Nfft  : Number of frequency bins.
%   gamma : Threshold on the STFT for reassignment,
%           default is 1E-6.
%
% OUTPUTS:
%   STFT  : The short-time Fourier transform.
%   SST   : Structure with vertical SST transforms for orders one and two.
%           They are respectively given by
%           SST.d1 and SST.d2
%   omega : Structure with IF estimations for orders one and two.
%           They are respectively given by
%           omega.d1 and omega.d2
%
% REFERENCES:
% [1] Behera, R., Meignen, S., & Oberlin, T. (2015). Theoretical Analysis
% of the Second-order Synchrosqueezing Transform. To appear in ACHA

if nargin == 3
    gamma = 1E-6;
end
 
 s = s(:);
 N = length(s);          
 
 ft   = 1:Nfft;
 bt   = 1:N;
  
 [g, l] = gauss_win(N, sigma);
 
 % Window definition
  
  n   = (0:2*l)'-l;
  t0  = n/N;
  t0  = t0(:);
  a   = pi/sigma^2;
  gp  = -2*a*t0.*g; 
  gpp = (-2*a+4*a^2*t0.^2).*g; % g''

 % Initialization
 STFT  = zeros(Nfft,N);
 Y  = zeros(Nfft,N);
 SST = struct('d1', Y, 'd2', Y);
 omega = struct('d1', Y, 'd2', Y);
 
 tau    = zeros(Nfft,N);
 phipp  = zeros(Nfft,N);
             
 %% Computes STFT and reassignment operators

 for b=1:N
 	% STFT, window g  
 	time_inst = -min([l,b-1]):min([l,N-b]);
    tmp = fft(s(bt(b)+time_inst).*g(l+time_inst+1),Nfft);
 	vg  = tmp(ft);
     
 	% STFT, window xg           
 	tmp = fft(s(bt(b)+time_inst).*(time_inst)'/N.*g(l+time_inst+1),Nfft);
 	vxg = tmp(ft);
       
    % operator Lx (dtau)
	
    tau(:,b) = vxg./vg;
 	
    % STFT, window gp
 	tmp = fft(s(bt(b)+time_inst).*gp(l+time_inst+1),Nfft);
 	vgp = tmp(ft);
    
    % operator omega
    omega.d1(:,b) = N/Nfft*(ft-1)'-real(vgp/2/1i/pi./vg);    
 	
    
    % STFT, window gpp
 	tmp  = fft(s(bt(b)+time_inst).*gpp(l+time_inst+1),Nfft);
 	vgpp = tmp(ft);
       
    %STFT, windox xgp
 	tmp  = fft(s(bt(b)+time_inst).*(time_inst)'/N.*gp(l+time_inst+1),Nfft);
 	vxgp = tmp(ft);
    
       
 	%computation of the two different omega 
        
    phipp(:,b) = 1/2/1i/pi*(vgpp.*vg-vgp.^2)./(vxg.*vgp-vxgp.*vg);
       
    %new omega2
    omega.d2(:,b) = omega.d1(:,b) - real(phipp(:,b)).*real(tau(:,b))...
                              + imag(phipp(:,b)).*imag(tau(:,b)); 

	% Storing STFT       
    STFT(:,b) = vg.*exp(2*1i*pi*(ft-1)'*min(l,b-1)/Nfft);%renormalized so that it fits with recmodes
 end
  
 %% reassignment step
 for b=1:N
    for eta=1:Nfft
        if abs(STFT(eta,b))> gamma
           k = 1+round(Nfft/N*omega.d1(eta,b));
            if (k >= 1) && (k <= Nfft)
             % original reassignment
             SST.d1(k,b) = SST.d1(k,b) + STFT(eta,b);
            end
            %reassignment using new omega2
            k = 1+round(Nfft/N*omega.d2(eta,b));
            if k>=1 && k<=Nfft
                % second-order Vertical reassignment: VSST
                SST.d2(k,b) = SST.d2(k,b) + STFT(eta,b);
            end 
        end
    end
 end
end