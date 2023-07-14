clear years

direc = '/home/disk/eos8/d.grosvenor/MODIS_L3_data/C6/';
%years = {'y2011_a','y2012_a','y2013_a','y2014_a'};
%years = {'y2011_a','y2012_a','y2013_a','y2014_a'};

years2=[2000:2014];
for i=1:length(years2)
    years{i}=['y' num2str(years2(i)) ''];
end

for iday=1:366
    for iyear=1:length(years)
        day = num2str(iday,'%03d');
        direc2 = [direc years{iyear} '/' day];
        files=dir(direc2);
        if length(files)==3 %dir returns . and .. as the first two
            cmd = ['!mv ' direc2 '/' files(3).name ' ' direc years{iyear} '/'];
            eval(cmd);
            cmd2 = ['!rmdir ' direc2];
            eval(cmd2);
        end
        
    end
    
end