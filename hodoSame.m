%plot hodograph U positive = moving towards the east and V positive = moving north - use readsound before




maxh=15000; %max height of hodograph in metres

heights=pr(2).p(:,1);

for i=1:length(pr)
    clear u ss
%sx=31; %spacing in u & v vectors
figure;
it=find(pr(i).p(:,1)>maxh);

%sx=round(it(1)/15);

for k=1:length(heights)
    sss=find( (pr(i).p(:,1)>=heights(k) |  abs(pr(i).p(:,1)-heights(k))<0.01 ) & pr(i).p(:,1)<maxh);
    if length(sss)>=1;ss(k)=sss(1);end
end

%u=-pr(i).p(1:sx:it(1),6);
%v=-pr(i).p(1:sx:it(1),5);
%h=pr(i).p(1:sx:it(1),1);

u=-pr(i).p(ss,6);
v=-pr(i).p(ss,5);
h=pr(i).p(ss,1);

[th,r]=cart2pol(u,v);

polarDan(th,r);
hold on;
plot(u,v,'ko');
x=[-25:0.1:25];
line(x,0);
line(0,x);

for j=1:length(u)
     text(u(j)+0.5,v(j)+0.1,num2str(h(j)/1000,'%.1f'),'fontweight','bold');
     fprintf(1,'%d %f\n',i,pr(i).p(sx*(j-1)+1,1));
end
axis equal; 

end
 