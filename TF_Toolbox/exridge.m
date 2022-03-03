function [Cs, Es] = exridge(TFR,jump,varargin)
% EXRIDGE extracts the ridges of a multicomponent signal
% [Cs, Es] = EXRIDGE(TFR,jump)
% [Cs, Es] = EXRIDGE(TFR,jump,nr,clwin)
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


