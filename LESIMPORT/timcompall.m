exdir='field\24.02.2015alt\timeseries\';
%ovride=1;
%nplots=2;
for i=4:56
    %con=i;
    [ha hb]=timcomp(SerDan,direcDan,textdataDan,3,i)


%gcf=h1;
%set(gcf,'paperpositionmode','auto');
exname=strcat('c:\matlabR12\work\',exdir,textdataDan(1).text,'-','col-',int2str(i),'-A.jpg');
print(ha,'-djpeg','-r150',exname);
close(ha);

exname=strcat('c:\matlabR12\work\',exdir,textdataDan(1).text,'-','col-',int2str(i),'-B.jpg');
print(hb,'-djpeg','-r150',exname);
close(hb);

end
%ovride=0;