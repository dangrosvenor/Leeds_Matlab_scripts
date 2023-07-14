outfile='c:\cygwin\home\user\lesexe\les2.3\tracer06-03-04\debugtwomey\SatDist';

fid=fopen(outfile,'w');

fprintf(fid,'&SatDist\n');

for i=1:length(snew)
    fprintf(fid,'SuperSats(%d)=%e,\n',i,snew(i));
end

for i=1:length(NnewT)
    fprintf(fid,'Naerosols(%d)=%e,\n',i,NnewT(i));
end

for i=1:length(NnewT)
    fprintf(fid,'MaeroAV(%d)=%e,\n',i,mmav2(i));
end

fprintf(fid,'&END');

fclose(fid);