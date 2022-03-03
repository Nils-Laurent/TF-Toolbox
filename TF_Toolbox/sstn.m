function  [STFT,SST,SST2,SST3,SST4,omega,omega2,omega3,omega4] = sstn(s,sigma,Nfft,gamma)

%% sstn_new : computes the STFT of a signal and different versions of synchrosqueezing
%
% INPUTS:   
%   s: real or complex signal
%   sigma: the parameter sigma in the definition of the Gaussian window  
%   Nfft: number of frequency bins
%   gamma: threshold on the STFT for reassignment
% OUTPUTS:   
%   STFT: the short-time Fourier transform
%   SST: standard synchrosqueezing
%   SST2: vertical second-order synchrosqueezing
%   SST3: vertical third-order synchrosqueezing
%   SST4: vertical fourth-order synchrosqueezing
%   omega: instantaneous frequency (vertical reassignment operator)
%   omega2: second-order instantaneous frequency
%   omega3: third-order instantaneous frequency
%   omega4: third-order instantaneous frequency
% REFERENCES
%[1] D.-H. Pham and S. Meignen, “High-order synchrosqueezing transform for 
%multicomponent signals analysis - with an application to gravitational-wave signal,”
%IEEE Transac tions on Signal Processing, vol. 65, pp. 3168–3178, June 2017.

s = s(:);

N = length(s);
ft   = (0:Nfft-1)*N/Nfft;
prec = 10^(-3);
L =  sigma*N;
l = floor(L*sqrt(-log(prec)/pi));
g = amgauss(2*l+1,l+1,L);

% Window definition

n   = (0:2*l)'-l;
t0  = n/N;
t0  = t0(:);
a   = pi/sigma^2; 
gp  = -2*a*t0.*g; 

% Initialization
STFT = zeros(Nfft,N);
SST = zeros(Nfft,N);
SST2 = zeros(Nfft,N);
SST3 = zeros(Nfft,N);
SST4 = zeros(Nfft,N);
omega = zeros(Nfft,N);
tau2 = zeros(Nfft,N);
tau3 = zeros(Nfft,N);
tau4 = zeros(Nfft,N);
omega2 = zeros(Nfft,N);
omega3 = zeros(Nfft,N);
omega4 = zeros(Nfft,N);
phi22p = zeros(Nfft,N);
phi23p = zeros(Nfft,N);
phi33p = zeros(Nfft,N);
phi24p = zeros(Nfft,N);
phi34p = zeros(Nfft,N);
phi44p = zeros(Nfft,N);

vg = zeros(Nfft,7);
vgp = zeros(Nfft,5);
Y = zeros(Nfft,4,4);
             
 %% Computes STFT and reassignment operators
        
 for b=1:N
    time_inst = -min([l,b-1]):min([l,N-b]);
 	
    % STFT, window x^i*g  
    for i = 0:7        
     vg(:,i+1) = fft(s(b+time_inst).*((time_inst)'/N).^i.*g(l+time_inst+1),Nfft);
    end

    % STFT, window x^i*gp
    for i = 0:5
     vgp(:,i+1) = fft(s(b+time_inst).*((time_inst)'/N).^i.*gp(l+time_inst+1),Nfft);
    end

    %% second-order operator tau
    tau2(:,b) = vg(:,2)./vg(:,1);
    % third order operator tau
    tau3(:,b) = vg(:,3)./vg(:,1);
    % four order operator tau
    tau4(:,b) = vg(:,4)./vg(:,1);
        
     %% Y expressions
    for i = 1:7
        for j = 1:7
            if i>=j
                Y(:,i,j) = vg(:,1).*vg(:,i+1) - vg(:,j).*vg(:,i-j+2);
            end
        end
    end  
    
    %% W expressions
    W2 = -1/2/1i/pi*(vg(:,1).^2+vg(:,1).*vgp(:,2)-vg(:,2).*vgp(:,1));
    W3 = -1/2/1i/pi*(2*vg(:,1).*vg(:,2)+vg(:,1).*vgp(:,3)-vg(:,3).*vgp(:,1));
    W4 = -1/2/1i/pi*(2*vg(:,1).*vg(:,3)+2*vg(:,2).^2+vg(:,1).*vgp(:,4) - vg(:,4).*vgp(:,1)+vg(:,2).*vgp(:,3) - vg(:,3).*vgp(:,2));
    
    %% operator omega
    omega(:,b) = ft'-real(vgp(:,1)/2/1i/pi./vg(:,1));     
    
    %% operator hat p: estimations of frequency modulation  
    %SST2 
    phi22p(:,b) = W2./Y(:,2,2);
    omega2(:,b) = omega(:,b) - real(phi22p(:,b).*tau2(:,b));
        
    %SST3
    phi33p(:,b) = (W3.*Y(:,2,2)-W2.*Y(:,3,3))./(Y(:,4,3).*Y(:,2,2)-Y(:,3,2).*Y(:,3,3));   
    phi23p(:,b) = (W2./Y(:,2,2) - phi33p(:,b).*Y(:,3,2)./Y(:,2,2));
    omega3(:,b) = omega(:,b)-real(phi23p(:,b).*tau2(:,b))-real(phi33p(:,b).*tau3(:,b));
       
    %SST4      
    phi44p(:,b) =((Y(:,4,3).*Y(:,2,2)-Y(:,3,2).*Y(:,3,3)).*W4-(W3.*Y(:,2,2)-W2.*Y(:,3,3)).*(Y(:,5,4)+Y(:,5,3)-Y(:,5,2))+(W3.*Y(:,3,2)-W2.*Y(:,4,3)).*(Y(:,4,4)+Y(:,4,3)-Y(:,4,2)))...
       ./((Y(:,4,3).*Y(:,2,2)-Y(:,3,2).*Y(:,3,3)).*(Y(:,6,4)+Y(:,6,3)-Y(:,6,2))-(Y(:,5,3).*Y(:,2,2)-Y(:,4,2).*Y(:,3,3)).*(Y(:,5,4)+Y(:,5,3)-Y(:,5,2))+(Y(:,5,3).*Y(:,3,2)-Y(:,4,2).*Y(:,4,3)).*(Y(:,4,4)+Y(:,4,3)-Y(:,4,2)));
    phi34p(:,b) = (W3.*Y(:,2,2)-W2.*Y(:,3,3))./(Y(:,4,3).*Y(:,2,2)-Y(:,3,2).*Y(:,3,3))...
       -phi44p(:,b).*(Y(:,5,3).*Y(:,2,2)-Y(:,4,2).*Y(:,3,3))./(Y(:,4,3).*Y(:,2,2)-Y(:,3,2).*Y(:,3,3));
    phi24p(:,b) = W2./Y(:,2,2) - phi34p(:,b).*Y(:,3,2)./Y(:,2,2) - phi44p(:,b).*Y(:,4,2)./Y(:,2,2);
   
    omega4(:,b) = omega(:,b) - real(phi24p(:,b).*tau2(:,b))- real(phi34p(:,b).*tau3(:,b))- real(phi44p(:,b).*tau4(:,b));

	% Storing STFT    
    STFT(:,b) = vg(:,1).*(exp(2*1i*pi*(0:Nfft-1)'*min(l,b-1)/Nfft));   
 end
 
 %% reassignment step
for b=1:N
    for eta=1:Nfft
        if abs(STFT(eta,b))> gamma
         k = 1+round(Nfft/N*omega(eta,b));
         if (k >= 1) && (k <= Nfft)
          % original reassignment
          SST(k,b) = SST(k,b) + STFT(eta,b);
         end
         %reassignment using new omega2
         k = 1+round(Nfft/N*omega2(eta,b));
         if k>=1 && k<=Nfft
          % second-order reassignment: SST2
          SST2(k,b) = SST2(k,b) + STFT(eta,b);
         end
         %reassignment using new omega3
         k = 1+round(Nfft/N*omega3(eta,b));
         if k>=1 && k<=Nfft
          % third-order reassignment: SST3
          SST3(k,b) = SST3(k,b) + STFT(eta,b);
         end
         %reassignment using new omega4
         k = 1+round(Nfft/N*omega4(eta,b));
         if k>=1 && k<=Nfft
          % fourth-order reassignment: SST4
          SST4(k,b) = SST4(k,b) + STFT(eta,b);
         end
        end
    end
 end
end

