%plot hodograph U positive = moving towards the east and V positive = moving north - use readsound before




maxh=20000; %max height of hodograph in metres

    clear u
%sx=31; %spacing in u & v vectors
figure;

hs=[620 1000 2500:1500:20e3];

H=pr(1).p(:,1);
%H=GridDan(1).Z+620;

%it=find(pr(i).p(:,1)>maxh);
for i=1:length(hs)
	[mindiff inds(i)]=min(abs(H-hs(i)));
end
%it=find(GridDan(1).Z+620>maxh);


%sx=round(it(1)/13);
%inds=1:sx:it(1);

%u=-pr(i).p(inds,6);
%v=-pr(i).p(inds,5);


u=velsU(inds); %sin
v=vels(inds); %cos

%u=GridDan(1).UBAR(inds);
%v=GridDan(1).VBAR(inds);



[th,r]=cart2pol(u,v);

polarDan(th,r);
hold on;
plot(u,v,'ko');
x=[-25:0.1:25];
line(x,0);
line(0,x);

for j=1:length(u)
     text(u(j)+0.5,v(j)+0.1,num2str(H(inds(j))/1000,'%.1f'),'fontweight','bold');
%     fprintf(1,'%d %f\n',i,pr(i).p(sx*(j-1)+1,1));
      fprintf(1,'%d %f\n',i,H(inds(j))/1000,1);

end
axis equal; 

 