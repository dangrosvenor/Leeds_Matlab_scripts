%loads radar stats and runs getstats to get rad strucure for radar and les data. hrstart is the start time  required for radar data
function [rad,iecho,date,day,month,year,hrs,mins,hrstart2,time,ia,ib,a,b,hrstartles]=loadRadStats(itype,Surf,hrstartles,timstart,timend)

time=[-300:300:15*3600];


switch itype
    case 1
		load 'C:\matlabR12\work\bauru\casestudy\radarSTATS\24.02\echo_3-end';
		load 'C:\matlabR12\work\bauru\casestudy\radarSTATS\24.02\timesEcho';
    case 2
        load 'C:\matlabR12\work\bauru\casestudy\radarSTATS\24.02\ppi_3-end';
		load 'C:\matlabR12\work\bauru\casestudy\radarSTATS\24.02\timesPPI_3-end';
    case 3
        load 'C:\matlabR12\work\bauru\casestudy\radarSTATS\24.02\cappi_3-end';
		load 'C:\matlabR12\work\bauru\casestudy\radarSTATS\24.02\timesCappi';
    case 4
        load 'C:\matlabR12\work\bauru\casestudy\radarSTATS\24.02\ppi_3-end';
		load 'C:\matlabR12\work\bauru\casestudy\radarSTATS\24.02\timesPPI_3-end';
    case 5
        load 'C:\matlabR12\work\bauru\casestudy\radarSTATS\24.02\cdmaxecho'; %load max echo radar data
                
end

if timstart==99
    timstart=8.9;
    timend=22.9;
end
    

[temp ia]=min(abs(time-(timstart-hrstartles)*3600));
[temp ib]=min(abs(time-(timend-hrstartles)*3600));

a=find(hrs>=timstart); %timstart=time for start of radar data required
a=a(1); %index for start time required
b=find(hrs>=timend); %timstart=time for start of radar data required
b=b(1); %index for start time required

iecho=itype;
rad=getstats(itype,Surf,np,date,day,month,year,hrs,mins,ia,ib,a,b);
hrstart2=timstart;
%hrstartles2=hrstartles;