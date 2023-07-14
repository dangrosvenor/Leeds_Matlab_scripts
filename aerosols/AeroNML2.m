outfile='c:\cygwin\home\user\lesexe\les2.3\tracer06-03-04\debugtwomey\MassDist';

fid=fopen(outfile,'w');

fprintf(fid,'&MassDist\n');

M2=fliplr(M);
for i=1:length(M)
    fprintf(fid,'MassAero(%d)=%e,\n',i,M2(i));
end

N2=fliplr(N);
for i=1:length(N)
    fprintf(fid,'Naerosols(%d)=%e,\n',i,N2(i));
end

Mav2=fliplr(Mav);
for i=1:length(Mav)
    fprintf(fid,'MaeroAV(%d)=%e,\n',i,Mav2(i));
end

fprintf(fid,'&END');

fclose(fid);