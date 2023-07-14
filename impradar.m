clear x a2 map c rhs time date

di1='c:/documents and settings/user/my documents/hibiscus/case_study/gerharddata/';
cd(di1);
d=dir;
[a b]=size(d);

for i=1:a
    k=findstr('top',d(i).name);
    [a2(i) b2]=size(k);    
end

[b ki]=find(a2==1); %only filenemes with 'top' in.
[c2 b2]=size(ki);

break;
for ii=1:b2
    date(ii).d=d(ki(ii)).name(7:13);
    day=str2num(date(ii).d(1:2));
    hrs(ii).t=str2num(d(ki(ii)).name(15:16));
    mins(ii).t=str2num(d(ki(ii)).name(17:18));
    
    
    fi=strcat(di1,d(ki(ii)).name);
    [im map]=imread(fi,'gif');
    

    
		if ii==1
			i=547;
			cnew=0;
			cold=0;
			
			values=[1.75:0.5:17.25]; %values of radar scale
			
			for j=1:32
                while cnew==cold  %until colour changes
                    cold=cnew;
                    i=i-1;   %move up image
                    cnew=im(i,586); %postion of colour scale on image
                end
                c(j)=cnew;  %store new colour
                cold=cnew;
			end %finds colour numbers of radar scale
			
			%c(33)=0; %black gridlines-put as NaN as don't know values
			%values(33)=0; %black	
            
        end

	for i=65:542  %just use square containing radar field
          for j=13:490
              rh=find(im(i,j)==c(:));
              aa=size(rh);
              if aa==0 | rh>30  %if height above 16.15km discard (not defined by radar)
                  rhs(ii).r(i-64,j-12)=NaN;
              else
                  
                 rhs(ii).r(i-64,j-12)=values(rh);
              end
          end
      end
 
       
end
save radarTOPS rhs time 
  
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

 
         
           


    

