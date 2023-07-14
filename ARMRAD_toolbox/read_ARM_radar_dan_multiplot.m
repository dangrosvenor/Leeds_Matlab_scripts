filePATH='z:/ACTIVE/PROCESSED/CURRENT/BOM-RADAR/20051116/';

lis=dir(filePATH);

%clear lis
%lis(3).name='cpol_hydroclass_20051116_0800.ascii';



%lis(3).name='filePATH

iplotandsave=0;

for iarm_radar=3:3%length(lis)
    fileNAME=lis(iarm_radar).name;
    if ~strcmp(fileNAME(end-3:end),'.zip')
    [xar,yar,zar,zh,hclass,slatr,slongr,adate,ahhmm,latitude,longitude,radClass]=...
        radarClassARMPaul([filePATH fileNAME]);
	end
    
    iarm_radar
    
    echo_tops;
    
    if iplotandsave==1
    
                    exdirA=['g:/ARM_radar/echo_tops/20051116/'];

                    bigcbar=0; %flag for one big colorbar instead of one for each subplot
					isamescale2=0;
					subplotting=0; %flag for whether to use subplots instead of new figures
                    lememm=0;
                    
					plotTimeHeightVap3;
                    jc=1;
                    exname=strcat(exdirA,fileNAME(26:29),'.emf');
					set(hf,'paperpositionmode','auto');
					%print(gcf,'-djpeg','-r350',exname);
                    
					print(hf,'-dmeta',exname);
					close(hf);
    end
                    
end

    
'done'