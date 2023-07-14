ccnflag=0;
figure('name','chaintrap24Mar');
nplots=3;
rowtot=4;
ccnchar=14;
maxav='avr';

rowno=1;
con=19;
maxav='max';
ytit='Max Total Rain Mixing Ratio (kg/m^2)';
trap22;

rowno=2;
con=18;
ytit='Average Total Ice Mixing Ratio (kg/m^2)';
maxav='avr'
trap22;

rowno=2;
con=7;
ytit='Time Average of Average Vertical Velocity (m/s)';
%trap22;

rowno=3;
con=6;
ytit='Average Vertical Mass Flux (kg/m^2/s)';
trap22;

rowno=4;
con=8;
ytit='Max Height of Threshold Ice Mixing Ratio (m)';
maxav='max';
%trap22;

rowno=4;
fnmax=15;
thresh=1e-5;
TRADD=3;
ytit='Max Height Reached by Tracer (m)';
heightice2;
%trap22;
for j=1:nplots
    dat=iceheights(j).SER;
    con=2;
    TrapIntDatN2;
end;
figlab='Averages of max heights of any ice';

ccnplN22;




axis(hs(1).h);
%axis(hs(1).h,[0 1.44e9 0.35 ans(4)]);
axis(hs(3).h);
%axis(hs(3).h,[0 1.44e9 0.085 ans(4)]);
axis(hs(4).h);
axis(hs(4).h,[0 1.44e9 7500 ans(4)]);
%axis(hs(5).h);
%axis(hs(5).h,[0 1.44e9 8000 ans(4)]);

axes('Position',[0 0 1 1],'Visible','off');
%text(.92,.22,'(e)');
%text(.92,.39,'(d)');
%text(.92,.56,'(c)');
%text(.92,.73,'(b)');
%text(.92,.9,'(a)');
