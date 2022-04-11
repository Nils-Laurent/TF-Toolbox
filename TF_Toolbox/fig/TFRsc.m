function TFRsc(varargin)

p = inputParser;

if nargin == 1
  p.addRequired(p, 'TFR');
  parse(p, varargin);
  TFR = p.Results.TFR;
  XV = 1:size(TFR, 2);
  YV = 1:size(TFR, 1);

  TFRsc(XV, YV, TFR, "", "");
end

p.addRequired(p, 'xVec');
p.addRequired(p, 'yVec');
p.addRequired(p, 'TFR');
p.addOptional(p, 'xunit', "time");
p.addOptional(p, 'yunit', "frequency");
p.addParameter(p, 'flipCBar', 1);
p.addParameter(p, 'axisFSZ', 18);
p.addParameter(p, 'labelFSZ', 18);
p.addParameter(p, 'squareSZ', 500);

parse(p, varargin);
r = p.Results;

if flipCBar
    cbar = flipud(gray);
else
    cbar = gray;
end

figure;
imagesc(r.xVec, r.Yvec, r.TFR);
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

