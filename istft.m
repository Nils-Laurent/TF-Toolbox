function [x] = istft(tfr,h,varargin)
%ISTFT inverse short time Fourier transform
%   [x] = istft(s,h)
%   tfr      : stft of x
%   h        : the filter
%
%   --  List of possible name-value pair argument
%   'cas'    : 1, no assumption on h except it does not vanish at 0
%              2, use a filter with unit energy
%              3, use a filter with unit mean
%              default is 1.
%   'len'    : length of the signal to be reconstructed.
%              default is size(s, 2).
%   'shift'  : parameter used to shift the downsampling
%              default is 0.
%
%   --  Output
%   x        : reconstruction from tfr
 
defaultCas = 1;
defaultLen = size(tfr, 2);
defaultShift = 0;

p = inputParser;
addRequired(p,'tfr');
addRequired(p,'h');
addParameter(p,'cas',defaultCas);
addParameter(p,'len',defaultLen);
addParameter(p,'shift',defaultShift);
parse(p,tfr,h,varargin{:});
cas = p.Results.cas;
Nsig = p.Results.len;
shift = p.Results.shift;

 [N,xrow_down] = size(tfr);
 downsamp = floor(Nsig/xrow_down);
 
 if (cas == 1)
  %case without periodizing
  x = zeros(Nsig,1);   
  Lh = (length(h)-1)/2;
  for icol=1:xrow_down,
   for q = 1:downsamp,
    ind = (icol-1)*downsamp+q; %mR+q
    if (ind <= Nsig)
     x(ind) = 1/h(Lh+q-shift)*mean(tfr(:,icol).*exp(2*pi*1i*(q-1-shift)*(0:N-1)'/N));
    end
   end
  end  
 end
 
 if (cas == 2)
  %case with periodization and with unit energy filters   
  %the filter h is with unit energy on its support
  Lh = (length(h)-1)/2; %h is supposed to be with odd length
  x = zeros(Nsig,1);
  for i = 1:Nsig   
   ind  = 1+(floor(((i-1)-Lh-shift)/downsamp):ceil(((i-1)+Lh+shift)/downsamp));
   ind = ind(((i-1)-downsamp*(ind-1)-shift >= -Lh)&((i-1)-downsamp*(ind-1)-shift <= Lh));
   %we renormalize the filter so that it is with unit energy (on the
   %samples corresponding to the downsampling)
   if (i > Lh)&&(i <= Nsig -Lh)
    x(i) = mean((tfr(:,ind).*exp(2*1i*pi*(0:N-1)'*((i-1)-downsamp*(ind-1)-shift)/N))*...
           h(Lh+1+(i-1)-downsamp*(ind-1)-shift))/norm(h(Lh+1+(i-1)-downsamp*(ind-1)-shift))^2;
   else
    %we renormalize the filter so that it is with unit energy (on the
    %samples corresponding to the downsampling)
    %To work well downsamp has to divide Nsig (cf paper)
    x(i) = mean((tfr(:,1+rem((ind-1)+xrow_down,xrow_down)).*exp(2*1i*pi*(0:N-1)'*((i-1)-downsamp*(ind-1)-shift)/N))...
                *h(Lh+1+(i-1)-downsamp*(ind-1)-shift))/norm(h(Lh+1+(i-1)-downsamp*(ind-1)-shift))^2;     
   end 
  end
 end
 
 if (cas == 3)
  %case with periodization and with unit l1 norm filters   
  %the filter h is with unit l1 norm
  
  Lh = (length(h)-1)/2;
  x = zeros(Nsig,1);
  
  for i = 1:Nsig   
   ind  = 1+(floor(((i-1)-Lh-shift)/downsamp):ceil(((i-1)+Lh+shift)/downsamp));
   ind = ind(((i-1)-downsamp*(ind-1)-shift >= -Lh)&((i-1)-downsamp*(ind-1)-shift <= Lh));
   if (i > Lh)&&(i <= Nsig-Lh)
    x(i) = mean((tfr(:,ind).*exp(2*1i*pi*(0:N-1)'*((i-1)-downsamp*(ind-1)-shift)/N))*...
           ones(length(Lh+1+(i-1)-downsamp*(ind-1)),1))/sum(h(Lh+1+(i-1)-downsamp*(ind-1)-shift));
   else
    x(i) = mean((tfr(:,1+rem((ind-1)+xrow_down,xrow_down)).*exp(2*1i*pi*(0:N-1)'*((i-1)-downsamp*(ind-1)-shift)/N))*...
           ones(length(Lh+1+(i-1)-downsamp*(ind-1)),1))/sum(h(Lh+1+(i-1)-downsamp*(ind-1)-shift));     
   end 
  end 
 end
end
