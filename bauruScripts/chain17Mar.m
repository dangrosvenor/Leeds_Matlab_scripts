figure('name',strcat('Column:',cons,' v.time'),'Position',posit);
nplots=6;

rowno=1;
tits;
ytit='Total Liquid Water Mixing Ratio (kg/kg)';
con=17;
timcomp2;

rowno=2;
tits2;
con=23
ytit='Max Rain Mixing Ratio (kg/kg)';
timcomp2;


rowno=3;
con=18;
ytit='';
timcomp2;


rowno=3;
con=18;
ytit='Total Ice Mixing Ratio (kg/kg)';
timcomp2;


rowno=4;
con=25;
ytit='Max Graupel Mixing Ratio (kg/kg)';
timcomp2;

rowno=5;
con=13;
ytit='Total Precipitation Rate (mm/hour)';
timcomp2;


axes('Position',[0 0 1 1],'Visible','off');
text(.92,.22,'(e)');
text(.92,.39,'(d)');
text(.92,.56,'(c)');
text(.92,.73,'(b)');
text(.92,.90,'(a)');
