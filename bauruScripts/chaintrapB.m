maxflag=0;
totflag=0;
divfac=1;

scrsz=get(0,'ScreenSize');
posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.1];
hf=figure('position',posit);

nplots=12;
rowtot=3;

rowno=1;
clear cono;
%cono=[22 24 25]; %eg for sum of cols 22,24&25
cono=[22 24 25];

exdir='bauru\homogen\zurich\avs2'; %filepath to export to


%ytit='Average Autoconversion Process Rate (kg/m^2/s)';
%ytit='Average Liquid Water Mixing Ratio (kg/m^2)';
%ytit='Average Rain Mixing Ratio (kg/m^2)';
ytit='Average Total Ice+Snow+Graupel Mixing Ratio (kg/m^2)';
%ytit='Height average (between 0 and 4000m) of time averaged (between 1000 and 2000s) profile of updraught velocity (m/s)';
%ytit='Height average of time averaged (between 1000 and 2000s) profile of updraught velocity (m/s)';
%ytit='Time Integrated Total Ice+Snow+Graupel Mixing Ratio Source Rates (kg/m^2)';
%ytit='Total precipitation at ground over 4000s (mm)';


maxflag=0;
trap2b;  %for timeseries averages
%trap2bpro; %for profile averages

rowno=2;
cono=6;
ytit='Average Vertical Mass Flux (kg/m^2/s)';
trap2b;

rowno=3;
cono=17;
ytit='Maximum Height of Cloud (km)';
maxflag=1;
divfac=1000;
trap2b;

maxflag=0;
divfac=1;

% rowno=4;
% cono=17;
% ytit='Average Height of Cloud Top (km)';
% trap2b;

%set(hs(1).h,'ylim',[0.069 0.071]);
%set(hs(1).h,'ylim',[0.265 0.275]);

%axis(hs(1).h);
%axis(hs(1).h,[0 960 0.12 ans(4)]);
% axis(hs(3).h);
% axis(hs(3).h,[0 960 10.5 ans(4)]);
% axis(hs(4).h);
% axis(hs(4).h,[0 960 8 ans(4)]);
%axis(hs(5).h);
%axis(hs(5).h,[0 1.44e9 8000 ans(4)]);

axes('Position',[0 0 1 1],'Visible','off');
%text(.92,.22,'(e)');
%text(.92,.39,'(d)');
%text(.92,.56,'(c)');
%text(.92,.73,'(b)');
%text(.92,.9,'(a)');



exname=strcat('c:\matlab6p1\work\',exdir);   
set(hf,'paperpositionmode','auto');
print(hf,'-djpeg','-r150',exname);

