for idat=1:length(ilats)
    ilat=ilats(idat);
    ilon=ilons(idat);
    
    m_plot(LONS(ilon)+0.5,LATS(ilat)+0.5,'ko','markerfacecolor','k');
    
    meanNd_mockL3(ilat,ilon)
    
    
end