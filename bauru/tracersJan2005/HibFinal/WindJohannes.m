%exdir='c:/matlabr12/work/bauru/tracersJan2005/force+3_3th3qv/diags/WindDiags50-105.mat'

graph=3;

timesLT=times/3600 + 9;

switch graph

case 1

    con(1).c=[0:0.02:10]; 
	con(2).c=[0:0.04:10]; 
    
	for  i=1:length(diag)
		figure;
        dgfind=findhead('ALu_A ',dgstrDan(i).dg);
        area=squeeze(diag(i).dg(1:100,dgfind(1),50:end));
        a=find(area==0);
        area(a)=1;
        
        dgfind=findhead('ALu_W ',dgstrDan(i).dg);
        dat=squeeze(diag(i).dg(1:100,dgfind(1),50:end))./area;
        [a b]=contourf(times(50:end),Grid.Z(1:100),dat,con(i).c);
        colorbarf(a,b);
        
	end


case 2
    

	for  i=1:length(diag)
        
        con(1).c=[0:-0.25/10:-0.25]; 
	    con(2).c=[0:-0.6/10:-0.6]; 
    
		figure;
        dgfind=findhead('ALd_A ',dgstrDan(i).dg);
        area=squeeze(diag(i).dg(1:100,dgfind(1),50:end));
        a=find(area==0);
        area(a)=1;
        
        dgfind=findhead('ALd_W ',dgstrDan(i).dg);
        dat=squeeze(diag(i).dg(1:100,dgfind(1),50:end))./area;
%         maxdat=max(max(dat));
%         mindat=min(min(dat));
%         inter=(maxdat-mindat)/12;
%         if maxdat>mindat
%             inter=-inter;
%             temp=minval;
%             minval=maxval;
%             maxval=temp;
%         end
        
        [a b]=contourf(times(50:end),Grid.Z(1:100),dat,con(i).c);
        colorbarf(a,b);
        
    end
    
    
case 3
    
    
    con(1).c=[0:1:16]; 
	con(2).c=[0:2:30];
    
    maxHi=110;
    
    for i=1:length(diag)
        dat=squeeze(max(wind(i).W(1:maxHi,:,50:end),[],2));
        figure;
        [a b]=contourf(timesLT(50:end),Grid.Z(1:maxHi)/1000,dat,con(i).c);
        colorbarf(a,b);
        
        xlabel('Local Time');
        ylabel('Height (km)');
        title('Max Updraught (m/s)');
    
    end
    
    
        
end