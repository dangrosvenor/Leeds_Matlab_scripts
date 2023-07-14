function plotSou(pro,col,GridDan,time)

    labs(1).l='Corumba';
	labs(2).l='Campo Grande';
	labs(3).l='Sau Paulo';
	labs(4).l='Curitiba';
	labs(5).l='Bauru';
    


if length(col)==2
    
    tit=strcat('Station:',labs(col(2)).l);
    clear labs;
    
    xdat(1).x=pro(col(2)).p(:,col(1),1);
    ydat(1).y=pro(col(2)).p(:,1,1);
    
    xdat(2).x=pro(col(2)).p(:,col(1),2);
    ydat(2).y=pro(col(2)).p(:,1,2);
    
    xdat(3).x=pro(end).p(:,col(1),1);
    ydat(3).y=pro(end).p(:,1,1);
%     
%     xdat(4).x=pro(end).p(:,col(1),2);
%     ydat(4).y=pro(end).p(:,1,2);
%     
%     xdat(5).x=GridDan(1).OLQBAR(:,1);
%     ydat(5).y=GridDan(1).Z;
    
    labs(1).l='24th 12';
	labs(2).l='25th 00';
    labs(3).l='Bauru 24th 20:15';
%     labs(4).l='Bauru 26th 12';
%     labs(5).l='LEM after 5hrs';
    
    
    
    
else
    
    tit=strcat('Time',int2str(time));
    
    
if col==11
        for i=1:size(pro,2)
            
            th=(pro(i).p(:,3,time)+273.15).*(1000./pro(i).p(:,2,1)).^0.286;
            xdat(i).x=pro(i).p(:,3,time);
            ydat(i).y=th;
            
		end
    
else
	for i=1:size(pro,2)
        
        xdat(i).x=pro(i).p(:,col,time);
        ydat(i).y=pro(i).p(:,1,time);
        
	end
end    
    labs(1).l='Corumba';
	labs(2).l='Campo Grande';
	labs(3).l='Sau Paulo';
	labs(4).l='Curitiba';
	labs(5).l='Bauru';
	labs(6).l='Bauru2';
    
    
    
    
    
    
end



figure('name',tit);
plotXY(xdat,ydat,labs,20);

title(tit);
