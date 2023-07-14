function Tephiplot(Tdat,Pdat,RHdat,Zdat)

% Run variables
flag=0;
dp=-10; %pressure step


% Physical constants
A=2.53e11;
epsilon=0.622;
B=5.42e3;
k=0.286;
L=2.5e6;
cp=1005;

% work out theta
THdat=(Tdat).*(1e5./Pdat).^k;
% initiallise
pstart=Pdat(1);
tstart=Tdat(1);
pend=2000;
rhstart=RHdat(1); %relative humidity at ground (%)


T(1)=tstart;
p(1)=pstart;
thstart=T(1)*(1e5/p(1))^k;
iend=-fix( (pstart-pend)/dp )
esg=A*exp(-B/T(1)); %es at ground
wg=rhstart/100*esg*epsilon/p(1); %mixing ratio at ground=constant

psc=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Works out the theta,t of a parcel rising from ground
disp('Working out path of parcel');
for i=1:iend
    dwsdp=-epsilon*A*exp(-B/T(i))/p(i)/p(i);
    dwsdt=B*epsilon*A*exp(-B/T(i))/T(i)/T(i)/p(i);
    es=A*exp(-B/T(i)); %saturation vapour pressure at T
    e=p(i)*wg/epsilon;  %vapour pressure at T

    if e>es    %if saturation reached
    psc=psc+1;
        if psc==1
            psati=i;
        end
        dTdp=( (k/p(i)) - (L*dwsdp/T(i)/cp) ) / ( (1/T(i)) + (L*dwsdt/T(i)/cp) ); %psuedo adiabat
        T(i+1)=T(i)+dTdp*dp;
    else
        T(i+1)=thstart*(p(i)/1e5)^k;  %dry adiabat
    end
    p(i+1)=p(i)+dp;
end
th=T.*(1e5./p).^k;
disp('Done');
warning off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('name','Tephigram plot');
pcs=[1000 800 700 600 500 400 300 200 100 50]; %mb to display
ths=233:20:613;
thsfine=220:0.1:600;
ts=40:-20:-90;

minT=min(T);
maxT=max(T);
minTH=log10(thsfine(1));
maxTH=log10(thsfine(end));
Y=1/1.2*(maxT-minT); %ratio of page height to width for page transform

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots isobars and pseudo-adiabat
for i=1:length(pcs)
    thcs(i,:)=T.*(1e3./pcs(i)).^k; % defines the theta line for the adiabat
    hold on;
    
    size(thcs);siz=ans(2);

    y1d=Y*(log10(thcs(1,1))-minTH)/(maxTH-minTH); %in page co-ordinates
    y2d=Y*(log10(thcs(1,siz))-minTH)/(maxTH-minTH);

    thi=atan( (y2d-y1d) / ( T(siz)-T(1) ) );
    %thi=0;
    [xx,yy]=rotatedan2(T,log(th)/log(10),thi,Y,minT,minTH,maxTH); 
    %h=plot(xx,(yy),'r','linewidth',2); %pseudo-adiabat

    [x,y]=rotatedan2(T,log(thcs(i,:))/log(10),thi,Y,minT,minTH,maxTH);
    h=plot( x,(y) ); %plot isobars

    tol=0;
    ii=[];
    % finds position to put text on pressure level
    while(size(ii,2)<1)
        tol=tol+2;
        ii=find(abs(x-xx(3))<tol); %find x value equal to xx(3)
    end
    h=text(xx(3),y(ii(end)),num2str(pcs(i)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=plot(xx,(yy),'r','linewidth',2); %pseudo-adiabat


clear x y
size(ths);
% Now plot adiabats
for i=1:ans(2)
   [x,y]=rotatedan2(T,log(ths(i))/log(10),thi,Y,minT,minTH,maxTH);
   h=plot( x,(y) );
   h=text(x(fix(siz*15/20)),y(fix(siz*15/20)),num2str(ths(i)),'rotation',-45);
end


clear x y
size(thsfine);
siz=ans(2);
size(ts);
% now plot isotherms
for i=1:ans(2)
   [x,y]=rotatedan2(273.13+ts(i),log(thsfine)/log(10),thi,Y,minT,minTH,maxTH);
   h=plot( x,(y) ); %plot isotherms
   h=text(x(fix(siz*6/20)),y(fix(siz*6/20)),num2str(ts(i)));

end


axis( [xx(end)/1.1 xx(1)*1.1 yy(1)*1.1 yy(end)*1.1 ] ); %zoom into area of interest
%set(gca,'XTicklabel','');
%set(gca,'YTicklabel',''); %remove axis tick lables
% now plot sounding
[x,y]=rotatedan2(Tdat,log(THdat)/log(10),thi,Y,minT,minTH,maxTH);
h=plot( x,(y),'k-x','linewidth',2); %plot sounding

[a b]=size(th);
[a c]=size(THdat);
sati=find(abs(Pdat-850)<2);

clear tdiff;
for i=100:b %start at low pressure to avoid getting highest pressure value as cross-over point
    [a2 ay]=min(abs(th(i)-THdat));    
    tdiff(i-99)=abs(T(i)-Tdat(ay));
%     if abs(tdiff(i)-0)<1e-10 
%         tdiff(i)=NaN;
%     end
end

[a b]=min(tdiff);
b=b+99;
fprintf(1,'Cross-over pressure = %g mb\n',p(b)/100);
[c d]=min(abs(Pdat.*100-p(b)));
fprintf(1,'Cross-over height = %g m',Zdat(d));
set(gcf,'position',[500 200 400 400])
set(gcf,'papersize',8.5.*[1 2.1])
