%process Marshak data

%reff = column 7
reffA = segA(7,:);
i0 = find(reffA==0);
reffA(i0)= NaN;

Xbins = [0:2:30];
[qh] = ndHistc_run([reffA(:)], Xbins');
mid_bins = 0.5 * (Xbins(1:end-1) + Xbins(2:end));
figure
plot(mid_bins,qh);

%reff = column 7
reffB = segB(3,:);
inan = find(reffB<=0);
reffB(inan)= NaN;

Xbins = [0:2:30];
[qh] = ndHistc_run([reffB(:)], Xbins');
mid_bins = 0.5 * (Xbins(1:end-1) + Xbins(2:end));
figure
plot(mid_bins,qh);

