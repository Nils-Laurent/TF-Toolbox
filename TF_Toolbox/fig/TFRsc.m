function TFRsc(varargin)

p = inputParser;

if nargin == 1
  addRequired(p, 'TFR');
  parse(p, varargin{:});
  TFR = p.Results.TFR;
  XV = 1:size(TFR, 2);
  YV = 1:size(TFR, 1);

  TFRsc(XV, YV, TFR, "", "");
  return;
end

addRequired(p, 'xVec');
addRequired(p, 'yVec');
addRequired(p, 'TFR');
addParameter(p, 'flipCBar', 1);
addParameter(p, 'axisFSZ', 20);
addParameter(p, 'labelFSZ', 20);
addParameter(p, 'squareSZ', 500);
addParameter(p, 'nfig', []);
addOptional(p, 'xunit', "time", @(x) 1 > 0);
addOptional(p, 'yunit', "frequency", @(x) 1 > 0);


if nargin > 4
  Index = find(strcmp(p.Parameters,varargin{4}), 1);
  if ~isempty(Index)
    TFRsc(varargin{1:3}, "time", "frequency", varargin{4:end});
    return;
  end
end

parse(p, varargin{:});
r = p.Results;

if sum(sum(abs(imag(r.TFR)))) > 0
  TFR = abs(r.TFR);
else
  TFR = r.TFR;
end

if r.flipCBar
    cbar = flipud(gray);
else
    cbar = gray;
end

if sum(size(r.nfig) > 0)
  figure(r.nfig);
else
  figure;
end
imagesc(r.xVec, r.yVec, TFR);
axis xy;
colormap(cbar);
ax = gca;
ax.XAxis.FontSize = r.axisFSZ;
ax.YAxis.FontSize = r.axisFSZ;
xlabel(r.xunit, 'interpreter', 'latex', 'FontSize', r.labelFSZ);
ylabel(r.yunit, 'interpreter', 'latex', 'FontSize', r.labelFSZ);
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, r.squareSZ, r.squareSZ]);

end

