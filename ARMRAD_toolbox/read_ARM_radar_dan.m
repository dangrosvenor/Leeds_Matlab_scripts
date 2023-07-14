filePATH='z:/ACTIVE/PROCESSED/CURRENT/BOM-RADAR/20051116/';
fileNAME='cpol_hydroclass_20051116_0400.ascii';

[xar,yar,zar,zh,hclass,slatr,slongr,adate,ahhmm,latitude,longitude,radClass]=...
        radarClassARMPaul([filePATH fileNAME]);
    
'done'