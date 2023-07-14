%correct slice data with values to be added on depending on radius

iload=0;
xc=333;
yc=305;

sfx=200/(533-333); %200km line situated at x=533 sf=1 - not necessary!

sfd(1)=200/(533-285); %0km at x=285 and 200km at x=533
sfd(2)=200/(550-295);
sfd(3)=200/(561-301);
sfd(4)=200/(488-269);
sfd(5)=100/(525-306);


xs(1)=119; %SECT_1_24
xs2(1)=506;
ys(1)=164;
ys2(1)=454;

xs(2)=333; %SECT_24
xs2(2)=333;
ys(2)=66;
ys2(2)=545;

xs(3)=280; %SECT_2_24
xs2(3)=518;
ys(3)=67;
ys2(3)=461;

xs(4)=340; %etc
xs2(4)=451;
ys(4)=66;
ys2(4)=461;

xs(5)=389;
xs2(5)=572;
ys(5)=70;
ys2(5)=274;


for i=1:length(xs)
    xs(i)=(xs(i)-xc);
    ys(i)=(yc-ys(i));
    xs2(i)=(xs2(i)-xc);
    ys2(i)=(yc-ys2(i));
end




rs=[15 40 90 115 160 200 240];
corrs=[0 4 8 8 9 11]; 
values=[10:1.5:56.5 999];

if iload==1
	'loading rad stats...'
		load 'C:\matlabR12\work\bauru\casestudy\radarSTATS\24.02\slices\sliceDat';    


end
'procressing...'
clear np4;

hdiff=hrs(2:end)-hrs(1:end-1);
fa=find(hdiff>2); %find where step up is large signifying the next slice - usually steps from 00hrs to 14hrs
fa(end+1)=length(aa);
fa(2:end+1)=fa(1:end);
fa(1)=0;



[x1 x2]=size(aa(1).i1);
for i=1:length(fa)-1
    np4(i).np(1:32,1:200)=0;
    x1=xs(i);
    x2=xs2(i);
    y1=ys(i);
    y2=ys2(i);
    theta=atan( (y1-y2) / (x2-x1) );
    th(i)=theta;
    for j=fa(i)+1:fa(i+1)
        if size(aa(j).i2,1)>0
            [ix iy]=find(aa(j).i2>0);
            %aa3(i).i2max(j-fa(i))=min(min(aa(j).i2(iz))); %only intersested in 10dBz echo and y direction (i2)
            for k=1:length(ix)
                d=(581-aa(j).i1(ix(k),iy(k)) )*sfd(i);
                yd=d*sin(theta);
                xd=d*cos(theta);
                xi=x1+xd;
                yi=y1-yd;
                r=sqrt(xi^2+yi^2);
                for kk=1:length(rs)-1
                    if ( r>rs(kk) & r<rs(kk+1) )
                        corr=corrs(kk);
                    end
                end
                dbz=values(ix(k))+corr;
                for kk=1:length(values)-1
                    if ( dbz>=values(kk) & dbz<values(kk+1) )
                        kval=kk;
                    end
                end
                np4(i).np(kval,j-fa(i))=np4(i).np(kval,j-fa(i))+1;
            end
        end    
    end                                
end