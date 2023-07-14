flt_graphs = {'Temp','Lat','Lon','Pressure','Altitude','Wind','Wind dir','Potemp'};

    itype='prof'; %(or timeseries)
    itimser = 'antjan06_flt';
    
for ifile=1:length(flt_graphs)
    flt_graph=flt_graphs{ifile};
    multisaveplot;
end