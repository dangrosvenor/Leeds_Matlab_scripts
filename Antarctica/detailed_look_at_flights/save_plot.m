

if isave_highlight==1    
    iremove=strfind(savename,':');
    savename(iremove)='_';
    
    iremove=strfind(savename,'\');
    savename(iremove)='';    
    
    iremove=strfind(savename,'/');
    savename(iremove)='';
    
    iremove=strfind(savename,'>');
    savename(iremove)='';
    
    iremove=strfind(savename,'<');
    savename(iremove)='';
    
    print(gcf,[savedir savename '.emf'],'-dmeta');    
%    print(gcf,[savedir savename '.eps'],'-depsc');    
    print(gcf,[savedir savename '.tiff'],'-dtiff');
    saveas(gcf,[savedir savename],'fig');
end

