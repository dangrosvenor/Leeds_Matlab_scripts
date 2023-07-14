fid=fopen('c:/temp/nitestNML','wt');

nitest=288;
sizex=192*3*500; %size of domain in x dir (m)
xtest=[-sizex/2:sizex/nitest:sizex/2]; %equally spaced points along x

fprintf(fid,'NITEST=%d,\n',length(xtest)); %note length of xtest = nitest +1 as have zero point in between
for i=1:length(xtest)
    fprintf(fid,'XTEST(%d)=%.1f,\n',i,xtest(i));
end

fclose(fid);

'finished writeXtestNML'