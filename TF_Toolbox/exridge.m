function [Cs, Es] = exridge(TFR,jump,varargin)
%EXRIDGE extracts the ridges of a multicomponent signal
% [Cs, Es] = EXRIDGE(TFR,jump)
% [Cs, Es] = EXRIDGE(TFR,jump,nr,clwin)
%
%   Extracts the ridge curve by maximising some energy.
%   The result is an approximation computed by a greedy algorithm.
%   The algorithm uses several random initializations and then
%   forward/backward search.
%
% INPUTS:
%   TFR     : time frequency representation
%   nr      : number of ridges
%             default is 1.
%
%   -- IF nr > 1, following argument has to be provided
%   clwin   : frequency clearing window.
%
% OUTPUTS:
%   Cs      : table of the nr ridges. Cs(:,j) is the j-th ridge location.

defaultNr = 1;
defaultClwin = -1;

p = inputParser;
addRequired(p,'tfr');
addRequired(p,'jump');
addOptional(p,'nr',defaultNr);
addOptional(p,'clwin',defaultClwin);
parse(p,TFR,jump,varargin{:});
nr    = p.Results.nr;
clwin = p.Results.clwin;

if nr > 1 && clwin < 0
    error("The clwin argument should be provided and has to be positive.");
end

[na,N] = size(TFR);

Cs = zeros(N, nr);
Es = zeros(nr, 1);

for j=1:(nr-1)
    [Cs(:,j), Es(j)] = exridge_mono(TFR,jump);

    % Remove this curve from the representation
    for b=1:N
        TFR(max(1,Cs(b,j)-clwin):min(na,Cs(b,j)+clwin),b)= 0;
    end
end
[Cs(:,nr), Es(nr)] = exridge_mono(TFR,jump);

Cs = Cs';

% sort the different ridges from HF to BF
[~, idx] = sort(sum(Cs,2),'ascend');
Cs = Cs(idx,:);

end

function [c,e] = exridge_mono(TFR,jump)
%EXRIDGE_MONO (internal usage) ridge detection
%   [c,e] = EXRIDGE_MONO(TFR,jump)
%
%   Extracts the ridge curve by maximising some energy.
%   The result is an approximation computed by a greedy algorithm.
%   The algorithm uses several random initializations and then
%   forward/backward search.
%
% INPUTS:
%   TFR    : time frequency representation
%   jump   : maximum difference for successive ridge points
%
% OUTPUTS:
%   c      : vector containing the indexes of the ridge
%   e      : energy of the returned ridge


%Et = log(abs(Tx)+eps^0.25);
Et     = abs(TFR).^2;
[na,N] = size(TFR);

% Parameters of the algorithm.
da = jump;
ng = min(60,floor(N/8)); % Number of random initializations.

ccur = zeros(1,N);
c = ccur;
e = -Inf;

cr_shift = @(v_idx, shift)(max(1,min(na,v_idx+shift)));
%     function [idx_shift] = internal_cr_shift(cr, )
%         idx_shift = max(1,idx-da):min(na,idx+da);
%     end

for k = floor(linspace(N/(ng+1),N-N/(ng+1),ng))
    [ecur,idx] = max(Et(:,k));
    ccur(k) = idx;
    Iq = max(1,idx-da):min(na,idx+da);
    [ecur,idx] = max(Et(Iq,k-1));
    ccur(k-1) = Iq(idx);

    cr = (ccur(k) - ccur(k-1));
%     idx = ccur(k) + cr;
    idx = cr_shift(ccur(k), cr);
    % forward step
    for b=k+1:N
     etmp = -Inf;
     a=max(1,idx-da):min(na,idx+da);

     [y, idx] = max(Et(a,b));
     idx = a(idx);
     ccur(b) = idx;
     etmp = y;
     ecur = ecur + etmp;

     cr = (ccur(b) - ccur(b-1));
%      idx = idx + cr;
     idx = cr_shift(idx, cr);
    end

    cr = (ccur(k+1) - ccur(k));
%     idx = ccur(k) - cr;
    idx = cr_shift(ccur(k), -cr);

    % backward step
    for b=k-1:-1:1
     etmp = -Inf;
     a=max(1,idx-da):min(na,idx+da);

     [y, idx] = max(Et(a,b));
     idx = a(idx);
     ccur(b) = idx;
     etmp = y;
     ecur = ecur + etmp;

     cr = (ccur(b+1) - ccur(b));
%      idx = idx - cr;
     idx = cr_shift(idx, -cr);
    end

    if ecur> e
        e = ecur;
        c = ccur;
    end
end

end

