%ppt radar dbZ stats 22/10/04


clear x a2 map c rhs time date dir vals zeros np aa hrs mins

%dire='c:/upload/radarims/ppi/24.02/';
dire='d:/';


if length(d)~=1896
    d=dir(dire);
end

[arad b]=size(d);

% for i=1:a
%     k=findstr('top',d(i).name);
%     [a2(i) b2]=size(k);    
% end

% [b ki]=find(a2==1); %only filenemes with 'top' in.
% [c2 b2]=size(ki);

dirout=strcat(dire,'stats');
%fid=fopen(dirout,'w');

%fprintf(fid,'%s %s %s %s %s %s %s %s %s %s %s %s ','date','time','13+km','13km','11.5km','10km','8.5km','7km','5.5km','4km','2.5km','1km');
%fprintf(fid,'%s %s %s %s %s %s %s %s %s %s\n','13+km','13km','11.5km','10km','8.5km','7km','5.5km','4km','2.5km','1km');

values=[8.5:1.5:55]; %values of radar scale (dbZ)

ii=0;
for ifile=1:arad
    
    if length(findstr('SECT_',d(ifile).name)>=1);
        ii=ii+1
    %if strcmp(d(ifile).name,'ppi_bru_240204_2100h.gif')==1
        
        fh=findstr('H.GIF',d(ifile).name);
        day(ii)=24;
        month(ii)=2;
        year(ii)=2004;
        hrs(ii)=str2num(d(ifile).name(fh-4:fh-3));
        %hrs=mod(hrs+3,24); %add 3 hrs to convert to UTC
        mins(ii)=str2num(d(ifile).name(fh-2:fh-1));
        
        
        fi=strcat(dire,d(ifile).name);
        [im map]=imread(fi,'gif');
    
        ix=760; %x pos of colorbar
        i=550; %y pos for bottom of cbar(first colour)
		cnew=999;
		cold=999;
		blackf=0;
		
        j=0;
		
		while cnew~=0 %until reach black
            while cnew==cold  %until colour changes
                cold=cnew;
                i=i-1;   %move up image
                cnew=im(i,ix); %postion of colour scale on image
            end
            j=j+1;
            c(j)=cnew;  %store new colour
            cold=cnew;
		end %finds colour numbers of radar scale
		
		%c(33)=0; %black gridlines-put as NaN as don't know values
		%values(33)=0; %black	
        
        
       
        
		im2=im(75:524,32:612); %y,x
        vals=zeros([size(im2,1) size(im2,2)]);
        zeroes=find(im2==0);  %points without black
        vals(zeroes)=NaN;
        %nzero=find(im2~=0);  %points without black
        for i=1:length(c)-1;
            a=find(im2==c(i));
            vals(a)=values(i);
        end
        
        %filterRadar; %filters out dodgy stuff in middle of PPIs
        
        for i=1:length(c)-1;
            a=find(vals==values(i));
            np(ii).np(i)=length(a);
            if length(a)>0
                aa(ii).a(i,1:np(ii).np(i))=a;
                aa(ii).i1(i,1:np(ii).np(i))=floor(a./size(im2,1))+1; %x
                aa(ii).i2(i,1:np(ii).np(i))=a-((aa(ii).i1(i,1:np(ii).np(i))'-1)*size(im2,1)); %y
            end 
        end
        
                    
% fprintf(fid,'%g/%g/%g %g:%g ',day,month,year,hrs,mins);                       
% for iout=1:size(dmin,2)
%     fprintf(fid,'%g ',dmin(iout)); %each pixel approx 1kmx1km
% end
% for iout=1:size(dmin,2)
%     fprintf(fid,'%g ',direc(iout));
% end
% fprintf(fid,'\n');

            
 
    end %if strcmp    
end %for ii=2:arad


%fclose(fid);
%save radarTOPS rhs time 
 


         
           


    

