for j=1:nplots
    clear dat;
    dat(:,1)=SerDan(j).SER(:,1);
    dat(:,2)=SerDan(j).SER(:,cono(1));
    size(cono);
    for ic=2:ans(2)
        dat(:,2)=dat(:,2)+SerDan(j).SER(:,cono(ic));
    end
    
%     if cono==99
%     dat(:,2)=SerDan(j).SER(:,22)+SerDan(j).SER(:,24)+SerDan(j).SER(:,25);    
%     else
%     dat(:,2)=SerDan(j).SER(:,cono);
%     end

    con=2;
    TrapIntDatN2;
end;
figlab=strcat('Averages of Time Series for column: ',int2str(con));
ccnplN2b;