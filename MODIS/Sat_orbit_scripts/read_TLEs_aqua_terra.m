%Read TLE archive files for AQUA and TERRA

filedir='C:\Users\Dan\Documents\logbook\UW\MODIS\Orbit prediction\TLE archive\';

filename_tle{1} = 'aqua/aqua_2002-2004_sat27424.txt';
filename_tle{2} = 'aqua/aqua_2005-2012_sat27424.txt';
filename_tle{3} = 'terra/terra_1999-2004_sat25994.txt';
filename_tle{4} = 'terra/terra_2005-2012_sat25994.txt';

clear tle01 tle02 satrec
for isat=1:2 %2 satellites
    iday=1;
    for ifile=1:2

        fid = fopen([filedir filename_tle{2*(isat-1)+ifile}],'rt');
        %the second files have a header of 5 lines - skip these
        if ifile==2
            for iskip=1:5
                tline=fgetl(fid);
            end
        end


        while 1
            tline = fgetl(fid);
            if ~ischar(tline) | strcmp(tline,'<End of file>')==1, break, end
            tle01{isat}{iday}=tline;
            
            tline = fgetl(fid);            
            if ~ischar(tline), break, end
            tle02{isat}{iday} = tline;
                       
            iday=iday+1;
        end
        fclose(fid);

    end
end

tle01_aqua = tle01{1};
tle02_aqua = tle02{1};

for itle=1:length(tle01_aqua)
    %converts the TLE into the format used by the orbit prediction
    %function
    satrec_aqua(itle) = twoline2rvMOD(tle01_aqua{itle},tle02_aqua{itle});
end

tle01_terra = tle01{2};
tle02_terra = tle02{2};

for itle=1:length(tle01_terra)
    %converts the TLE into the format used by the orbit prediction
    %function
    satrec_terra(itle) = twoline2rvMOD(tle01_terra{itle},tle02_terra{itle});
end



save_filename = 'aqua_terra_TLEs_MATLAB.mat';
save([filedir save_filename],'tle01_aqua','tle02_aqua','tle01_terra','tle02_terra','satrec_aqua','satrec_terra');


