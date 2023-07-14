%read in the night-flight leg data from Chris T

file_loc = '/home/disk/eos1/d.grosvenor/modis_work/VOCALS/';
filename = 'Night_legs_inclouddata_stats.txt';
filepath = [file_loc filename];

fid = fopen(filepath,'rt');

nf_dat = fscanf(fid,'%f\n',[5 Inf]);

%columns
%1 cloud leg
%2 day (counted in days after 30th Sep, 2008)
%3 Lat
%4 Lon
%5 seconds after midnight UTC

%quick conversion to local time is to add (longitude/360*24) to UTC


nf_time_matlab = datenum('30-Sep-2008') + nf_dat(2,:) + nf_dat(5,:)/3600/24; 

nf_sza = sun_pos(nf_time_matlab,nf_dat(3,:),nf_dat(4,:));

D = distlatlon(nf_dat(3,2:end),nf_dat(4,2:end),nf_dat(3,1:end-1),nf_dat(4,1:end-1));

figure('name','SZA','Position',[9 50 1200 300]);
plot(nf_time_matlab,nf_sza,'bo-');
xlabel('Date (UTC)');
ylabel('SZA');
d=datenum('18-Oct-2008') + [0:31];
grid
set(gca,'xtick',d);
set(gca,'xlim',[d(1) d(end)]);
datetickzoom('x','dd','keeplimits');




