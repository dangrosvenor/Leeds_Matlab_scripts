function admultiF(Tdat,Pdat,rhstart,Hdat);
%takes one AED sounding and takes out the text, saves as .txt file and plots tephigram

teph=1; %flag for tephigrams 1=yes 0=no
%fla=1; %1=AED 2=old Bauru style

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
% if fla==2
%     FileName='c:/documents and settings/user/my documents/cprogs/0302052040UTC.TXT';
% else
%     FileName='c:/documents and settings/user/my documents/HIBISCUS/baurufield/soundings/ascii/04_02_06_12AED.txt';
% end

% if fla==4
%     Tdat=dat(:,3);
%     Pdat=dat(:,2)
%     rhstart=dat(1,4); %relative humidity at ground (%)
%     Hdat=dat(:,1);
% elseif fla==3
%     Tdat=pr(ii).p(pind(1):end,3);
%     Pdat=pr(ii).p(pind(1):end,2);
%     rhstart=pr(ii).p(pind(1),4); %relative humidity at ground (%)
%     Hdat=pr(ii).p(pind(1):end,1);
% elseif fla==2
%     %txt=fscanf(fid,'%s',71);
%     dat=fscanf(fid,'%g',[9 inf]);
%     Tdat=dat(5,:);
%     Pdat=dat(4,:);
%     rhstart=dat(6,1); %relative humidity at ground (%)
% else
%     jmin=1;
%     clear dat;
%     i=0;
%     iflag=1;
%     
%     while iflag==1 & i<1000
%         i=i+1;
%         for j=jmin:9
%             c=fscanf(fid,'%g',[1 1]);
%             [a b]=size(c);
%             if a==1 %if is a float
%                 dat(j,i)=c;
%             else
%                 iflag=0; %else have reached missing data -exit
%                 break;
%             end
%         end
%         jmin=1;
%         flagread=1;
%         count=0;
%         while flagread==1 & iflag==1
%             c=fscanf(fid,'%g',[1 1]);
%             [a b]=size(c);
%             
%             if a==0 %if there is text at the end of the row
%                 count=count+1;
%                 c=fscanf(fid,'%c',[1 1]);
%                 if count>100 %if end of file is reached without encountering missing data will keep repeating
%                     iflag=0;
%                     break;
%                 end
%             else
%                 dat(1,i+1)=c;
%                 flagread=0;
%                 jmin=2; %otherwise value from next row is scanned so start from j=2
%             end
%         end
%                 
%     end
%     
%     
%     
%     dat(:,i)=[]; 
%     for j=1:i-1  %write to new file
%         for ii=1:9
%             fprintf(fiid,'%g ',dat(ii,j));
%             j
%             %dat(ii,j)
%         end
%         fprintf(fiid,'\n');
%     end
%     
%     outname
%     
%     
%     
%     Tdat=dat(5,:);
%     Pdat=dat(4,:);
%     rhstart=dat(6,1); %relative humidity at ground (%)
%     
% end %else

if teph==1 %big if

THdat=(Tdat+273.13).*(1e3./Pdat).^k;
pstart=Pdat(1)*100;
tstart=Tdat(1)+273.13;
pend=2000;


% pstart=100000;
% tstart=20+273;

T(1)=tstart;
p(1)=pstart;
thstart=T(1)*(1e5/p(1))^k;
iend=-fix( (pstart-pend)/dp )
esg=A*exp(-B/T(1)); %es at ground
wg=rhstart/100*esg*epsilon/p(1); %mixing ratio at ground=constant

psc=0;

for i=1:iend
    
    if i>100000
        'Breaking out'
        break;
    end
    
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
    %i
end

p(i+1)=p(i)+dp;

end

th=T.*(1e5./p).^k;

hold on;

figure;
pcs=[1000 800 700 600 500 400 300 200 100 50]; %mb
ths=273:20:413;
thsfine=260:0.1:600;
ts=30:-20:-90;

minT=min(T);
maxT=max(T);
minTH=log10(thsfine(1));
maxTH=log10(thsfine(end));
Y=1/1.2*(maxT-minT); %ratio of page height to width for page transform


size(pcs)
for i=1:ans(2)
    
thcs(i,:)=T.*(1e3./pcs(i)).^k;

hold on;
size(thcs);
siz=ans(2);

y1d=Y*(log10(thcs(1,1))-minTH)/(maxTH-minTH); %in page co-ordinates
y2d=Y*(log10(thcs(1,siz))-minTH)/(maxTH-minTH);

thi=atan( (y2d-y1d) / ( T(siz)-T(1) ) );
%thi=0;


[xx,yy]=rotatedan2(T,log(th)/log(10),thi,Y,minT,minTH,maxTH); 
h=plot(xx,(yy),'r'); %plot pseudo-adiabat
set(h,'linewidth',2);

[x,y]=rotatedan2(T,log(thcs(i,:))/log(10),thi,Y,minT,minTH,maxTH);
h=plot( x,(y) ); %plot isobars

tol=0;
ii=[];
while size(ii,2)<1
    tol=tol+2
    ii=find(abs(x-xx(3))<tol); %find x value equal to xx(3)
end
h=text(xx(3),y(ii(end)),num2str(pcs(i)));


end




clear x y
size(ths);
for i=1:ans(2)
   [x,y]=rotatedan2(T,log(ths(i))/log(10),thi,Y,minT,minTH,maxTH);
   h=plot( x,(y) );
   h=text(x(fix(siz*15/20)),y(fix(siz*15/20)),num2str(ths(i)),'rotation',-19);
end



clear x y
size(thsfine);
siz=ans(2);
size(ts);
for i=1:ans(2)
   [x,y]=rotatedan2(273.13+ts(i),log(thsfine)/log(10),thi,Y,minT,minTH,maxTH);
   h=plot( x,(y) ); %plot isotherms
   h=text(x(fix(siz*6/20)),y(fix(siz*6/20)),num2str(ts(i)));

end


axis( [xx(end)/1.1 xx(1)*1.1 yy(1)*1.1 yy(end)*1.1 ] ); %zoom into area of interest
set(gca,'XTicklabel','');
set(gca,'YTicklabel',''); %remove axis tick lables




[x,y]=rotatedan2(273.13+Tdat,log(THdat)/log(10),thi,Y,minT,minTH,maxTH);
h=plot( x,(y),'kx'); %plot sounding

[a b]=size(th);
[a c]=size(THdat);
sati=find(abs(Pdat-850)<2);

clear tdiff;
for i=100:b %start at low pressure to avoid getting highest pressure value as cross-over point
    [a2 ay]=min(abs(th(i)-THdat));    
    tdiff(i-99)=abs(T(i)-273-Tdat(ay));
%     if abs(tdiff(i)-0)<1e-10 
%         tdiff(i)=NaN;
%     end
end

[a b]=min(tdiff);
b=b+99;
fprintf(1,'Cross-over pressure = %g mb\n',p(b)/100);
[c d]=min(abs(Pdat.*100-p(b)));
fprintf(1,'Cross-over height = %g m',Hdat(d));

end %big if