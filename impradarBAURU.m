clear x a2 map c rhs time date dir

dire='c:/documents and settings/user/my documents/hibiscus/baurufield/radar/echotops/23-29.02/';
d=dir(dire);
[arad b]=size(d);

% for i=1:a
%     k=findstr('top',d(i).name);
%     [a2(i) b2]=size(k);    
% end

% [b ki]=find(a2==1); %only filenemes with 'top' in.
% [c2 b2]=size(ki);

dirout='c:/documents and settings/user/my documents/hibiscus/baurufield/radar/echotops/23-29.02/stats';
fid=fopen(dirout,'w');

fprintf(fid,'%s %s %s %s %s %s %s %s %s %s %s %s ','date','time','13+km','13km','11.5km','10km','8.5km','7km','5.5km','4km','2.5km','1km');
fprintf(fid,'%s %s %s %s %s %s %s %s %s %s\n','13+km','13km','11.5km','10km','8.5km','7km','5.5km','4km','2.5km','1km');

for ii=3:30  %arad
    ii
    date(ii).d=d(ii).name(1:6);
    day=str2num(date(ii).d(5:6));
    month=str2num(date(ii).d(3:4));
    year=str2num(date(ii).d(1:2));
    hrs=str2num(d(ii).name(7:8)) + 3; %minus 3 hrs to convert to UTC
    hrs=mod(hrs,24);
    mins=str2num(d(ii).name(10:11));
    
    
    fi=strcat(dire,d(ii).name);
    [im map]=imread(fi,'gif');
    

    
% 		if ii==1
% 			i=547;
% 			cnew=0;
% 			cold=0;
% 			
% 			values=[1.75:0.5:17.25]; %values of radar scale
% 			
% 			for j=1:32
%                 while cnew==cold  %until colour changes
%                     cold=cnew;
%                     i=i-1;   %move up image
%                     cnew=im(i,586); %postion of colour scale on image
%                 end
%                 c(j)=cnew;  %store new colour
%                 cold=cnew;
% 			end %finds colour numbers of radar scale
% 			
% 			%c(33)=0; %black gridlines-put as NaN as don't know values
% 			%values(33)=0; %black	
%             
%         end

	%for i=65:542  %just use square containing radar field
    hvals=11:-1:2; %colour lables for 13+,13,11.5,10,8.5,7,5.5,4,2.5,1km
    
    c1=242;
    c2=239; % co-ords for centre of radar image
    hmax=0;
    dmin(1:size(hvals,2))=999; 
    direc(1:size(hvals,2))=999;
    for j=1:478
              for k=1:size(hvals,2);
                  rh=find(im(:,j)==hvals(k)); %find points with value hvals(k)
                  aa=size(rh);
                  if aa>0 
                    for ih=1:aa;
                            dist=sqrt( (c2-j)^2 + (rh(ih)-c1)^2 ); %work out how far from centre
                            if dist<dmin(k)
                                %height=im(rh(ih),j);
                                dmin(k)=dist;
                                 
                                xd=(j-c2);
                                yd=(c1-rh(ih));
                                direc(k)=atan( yd/xd )*180/pi; %direction in degrees
                                
%                                 if k==4
%                                     j
%                                     rh(ih)
%                                     xd
%                                     yd
%                                  end
                                
                                if direc(k)>=0 %convert to degrees in 360 circle
                                    if xd<0 | yd<0
                                        direc(k)=direc(k)+180; %angle between 180 and 270
                                    end
                                else
                                    if yd>=0 & xd<=0
                                        direc(k)=180+direc(k); %angle between 90 and 180 (remember angle is negative)
                                    else
                                        direc(k)=360+direc(k); %angle between 270 and 360 (remember angle is negative)
                                    end
                                end
                                
                                direc(k)=mod(90-direc(k),360); %add 90 to make 0deg north
                              
          
                            end
                    end
                  end
              end
              
          end %for j=1:478
                    
fprintf(fid,'%g/%g/%g %g:%g ',day,month,year,hrs,mins);                       
for iout=1:size(dmin,2)
    fprintf(fid,'%g ',dmin(iout)); %each pixel approx 1kmx1km
end
for iout=1:size(dmin,2)
    fprintf(fid,'%g ',direc(iout));
end
fprintf(fid,'\n');

            
 
       
end
fclose(fid);
%save radarTOPS rhs time 
  
break;


[a b]=size(im);
starty=303;
startx=252;
R=239;
ndeg=500;
nr=300;

ddeg=(pi/180)*360/(ndeg-1);
dr=R/(nr-1);

r=0;
deg=0;


%          x=startx+round(r*cos(deg));
%          y=starty+r*round(sin(deg));
         

rs=repmat([0:dr:R],ndeg,1);
degs=repmat([0:ddeg:2*pi],nr,1)';
y=startx + round( rs.*cos(degs) );
x=starty + round( rs.*sin(degs) );


cols=im(1).im(sub2ind(size(im(1).im),x,y));

break;

 
         
           


    

