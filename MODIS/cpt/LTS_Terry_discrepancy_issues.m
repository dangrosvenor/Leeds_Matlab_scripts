%calculate some straight longitude transects for comparison to Terry's data
%rather than using the watervap case 115

Terry_textfile = '/home/disk/eos1/d.grosvenor/Terry_LTS_data.txt';


%find time indices for ASON (SON + August)
it=find(daynum_timeseries3_ERAInt==214 | daynum_timeseries3_ERAInt==245 | daynum_timeseries3_ERAInt==275 | daynum_timeseries3_ERAInt==306);
LTS = meanNoNan(ecLTS_mon_Man_ALL(it,:,:),1);
%Plat are the face values for lat (data is 121 x 240 and Plat is also 121
%long and runs from 90 to -90 with 1.5 deg resolution)

%find the face closest to -20
[minval,i20]=min(abs(Plat--20));  %i20 = 74
figure
plot(Plon,LTS(i20-1,:),'k--');  %Plat(i20-1) = -18
hold on
plot(Plon,LTS(i20,:));          %Plat(i20) = -19.5
plot(Plon,LTS(i20+1,:),'b--');  %Plat(i20+1) = -21.0

leg_str = {'-18 S','-19.5 S','-21.0 S'};
legend(leg_str,'location','northwest');
set(gca,'xlim',[-140 -70]);
xlabel('Long');
ylabel('LTS (K)');


fid=fopen(Terry_textfile,'wt');
fprintf(fid,'Lon LTS_18S LTS_19.5S LTS_21S\n');

ilons = find(Plon>-141 & Plon < -69);

for i=1:length(ilons)
    fprintf(fid,'%f %f %f %f\n',Plon(ilons(i)),LTS(i20-1,ilons(i)),LTS(i20,ilons(i)),LTS(i20+1,ilons(i)));
end

fclose(fid);






