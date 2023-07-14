flag=0;

dp=-10; %pressure step



clear p T th thcs r theta

A=2.53e11;
epsilon=0.622;
B=5.42e3;
k=0.286;
L=2.5e6;
cp=1005;


%read in sounding
FileName='c:/documents and settings/user/my documents/cprogs/0302052040UTC.TXT';
fid=fopen(FileName,'rb');
txt=fscanf(fid,'%s',71);
dat=fscanf(fid,'%g',[12,108]);


Tdat=dat(6,:);
Pdat=dat(5,:);
THdat=(Tdat+273.13).*(1e3./Pdat).^k;

pstart=Pdat(1)*100;
tstart=Tdat(1)+273.13;
pend=2000;
rhstart=dat(7,1); %relative humidity at ground (%)

T(1)=tstart;
p(1)=pstart;
thstart=T(1)*(1e5/p(1))^k;
iend=-fix( (pstart-pend)/dp )
esg=A*exp(-B/T(1)); %es at ground
wg=rhstart/100*esg*epsilon/p(1); %mixing ratio at ground=constant

for i=1:iend
dwsdp=-epsilon*A*exp(-B/T(i))/p(i)/p(i);
dwsdt=B*epsilon*A*exp(-B/T(i))/T(i)/T(i)/p(i);
es=A*exp(-B/T(i)); %saturation vapour pressure at T
e=p(i)*wg/epsilon;  %vapour pressure at T


psc=0;

if e>es    
    psc=psc+1;
    if psc==1
        psati=i;
    end
    dTdp=( (k/p(i)) - (L*dwsdp/T(i)/cp) ) / ( (1/T(i)) + (L*dwsdt/T(i)/cp) );
    T(i+1)=T(i)+dTdp*dp;
else
    i
    T(i+1)=thstart*(p(i)/1e5)^k;  
end

p(i+1)=p(i)+dp;

end

th=T.*(1e5./p).^k;
%figure;
%plot(T,log(th));
hold on;

figure;
pcs=[1000 800 600 400 300 200 100 68.3 40 20 10]; %mb
ths=260:40:600;
thsfine=260:0.1:600;
ts=30:-10:-90;



size(pcs)
for i=1:ans(2)
thcs(i,:)=T.*(1e3./pcs(i)).^k;
%plot(T,log(thcs(i,:)));
hold on;
size(thcs);
siz=ans(2);
thi=atan( (100*log(thcs(1,siz)) - 100*log(thcs(1,1))) / ( T(siz)-T(1) ) );
thi=0;
%thi= atan( log(thcs(1,siz)) / T(siz) ) - atan(log(thcs(1,1)) / T(1) );
[xx,yy]=rotatedan(T,100*log(th),thi);
h=plot(xx,(yy),'r'); %pseudo-adiabat


[x,y]=rotatedan(T,100*log(thcs(i,:)),thi);
%r(i,:)=( (T.^2)+(log(thcs(i,:)).^2) ).^0.5;
%theta(i,:)=atan( log(thcs(i,:))./T );
h=plot( x,(y) );
%rotate(h,[0 0 1],a*180/pi);    
h=text(x(fix(siz*19/20)),y(fix(siz*19/20)),num2str(pcs(i)));
%rotate(h,[0 0 1],a*180/pi);
%plot( r(i,:).*cos(theta(i,:)-thi) , log(r(i,:).*sin(theta(i,:)-thi)) );

end




clear x y
size(ths);
for i=1:ans(2)
   [x,y]=rotatedan(T,100*log(ths(i)),thi);
   h=plot( x,(y) );
   %rotate(h,[0 0 1],a*180/pi);
   h=text(x(fix(siz*8/20)),y(fix(siz*8/20)),num2str(ths(i)),'rotation',-19);
   %rotate(h,[0 0 1],a*180/pi);
end



clear x y
size(thsfine);
siz=ans(2);
size(ts);
for i=1:ans(2)
   [x,y]=rotatedan(273.13+ts(i),100*log(thsfine),thi);
   h=plot( x,(y) );
   %rotate(h,[0 0 1],a*10*180/pi);
   h=text(x(fix(siz*3/20)),y(fix(siz*3/20)),num2str(ts(i)));
   %rotate(h,[0 0 1],a*180/pi);
end




thi


%rotate(h,[0 0 1],a*180/pi) 
size(T);
%axis( [xx(ans(2)) xx(1) yy(1) yy(ans(2)) ] );
%set(gca,'XTicklabel','');
%set(gca,'YTicklabel','');




[x,y]=rotatedan(273.13+Tdat,100*log(THdat),thi);
h=plot( x,(y),'kx'); %plot sounding

[a b]=size(th);
[a c]=size(THdat);
sati=find(abs(Pdat-850)<2);

clear tdiff;
for i=100:b
    [a2 ay]=min(abs(th(i)-THdat));    
    tdiff(i-99)=abs(T(i)-273-Tdat(ay));
%     if abs(tdiff(i)-0)<1e-10 
%         tdiff(i)=NaN;
%     end
end

[a b]=min(tdiff);
b=b+99;
fprintf(1,'Cross-over pressure = %g\n',p(b)/100);
[c d]=min(abs(Pdat.*100-p(b)));
fprintf(1,'Cross-over height = %g',dat(4,d));


break;
rr=( (T.^2)+(log(th).^2) ).^0.5;
angle=atan( log(th)./T );
xx=rr.*cos(angle-thi);
yy=log(rr.*sin(angle-thi));
plot(xx,yy,'--','r');
size(T);
axis( [xx(ans(2)) xx(1) yy(1) yy(ans(2)) ] );


