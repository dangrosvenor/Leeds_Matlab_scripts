%prepare pr from LES output for writing sounding

pr(3).p(:,1)=Grid.Z;
 pr(3).p(:,2)=Grid.PREFN/100;
 T=tempLES(Grid);
 pr(3).p(:,3)=T-273.15;
 pr(3).p(:,10)=Grid.OLQBAR(:,1);
 pr(3).p(:,8)=Grid.VBAR(:,1);
 
 
