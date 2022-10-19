function write_figfiles(fname)

fname2 = string(fname);

fname2 = join(fname2, '');

saveas(gcf, fname2, 'epsc');
savefig(fname2);
end

