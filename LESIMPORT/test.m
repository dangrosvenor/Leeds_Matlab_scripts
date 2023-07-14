iso3d;
exname=strcat('c:\matlabR12\work\',exdir,direcDan(i).dir(end-2:end),textdataDan(i).text,'.jpg');
print(gcf,'-djpeg','-r150',exname);
close(gcf);