scrsz=get(0,'ScreenSize');
exist posit;
a=ans;
if a==0;
    posit=[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2.2];
end;

exist nstart;
a=ans;
if a==0;
    nstart=1;
    nend=6;
end

hf=figure('name',strcat('ChainBauru'),'Position',posit);

nplots=12;
rowtot=1;

rowno=1;
titsbauruFull;
ytit='Total Liquid Water MR ( kg/m ^2 )';
con=21;
timcomp2b;

rowno=2;
tits2;
con=23;
ytit='Total Rain MR ( kg/m ^2 )';
timcomp2b;

rowno=3;
tits2;
con=26;
ytit='Total Precipitation MR ( mm/hour )';
timcomp2b;



rowno=4;
con=22;
ytit='Total Ice MR ( kg/m ^2 )';
timcomp2b;


rowno=5;
con=25;
ytit='Total Graupel MR ( kg/m ^2 )';
timcomp2b;

rowno=6;
con=24;
ytit='Total Snow MR ( kg/m ^2 )';
timcomp2b;


axes('Position',[0 0 1 1],'Visible','off');
text(.92,.22,'(e)');
text(.92,.39,'(d)');
text(.92,.56,'(c)');
text(.92,.73,'(b)');
text(.92,.90,'(a)');
