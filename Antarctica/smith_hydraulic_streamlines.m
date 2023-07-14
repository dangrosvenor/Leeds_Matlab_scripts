H0 = 6e3;
H0 = 3e3;
H0 = 2.9e3;
%H0 = 3051;
%H0 = 2.15e3;

Y=zz(1).z'*1000;
X=d_line*1000;

% u = nc{'U'}(time,:,y_inds,x_inds);
%             v = nc{'V'}(time,:,y_inds,x_inds);
%             w = nc{'W'}(time,:,y_inds,x_inds);
%             z_temp=(nc{'PH'}(time,:,y_inds,x_inds) + nc{'PHB'}(time,:,y_inds,x_inds) )./9.81;
%             z = 0.5.*(z_temp(1:end-1,:,:)+z_temp(2:end,:,:));
%             terr_all = nc{'HGT'}(time,:);
%             terr = terr_all(y_inds,x_inds);
cont_data = nc{'T'}(time,:,y_inds,x_inds) + 300; %potemp


clear potemp u_section v_section

for i_line=1 % only need to do for the first profile of the cross-section
    %pdat(1).p(:,i_line-n_line(1)+1) = interp1 (z_slice(i_line,:)/1e3,v_slice(i_line,:),zz(1).z); %interpolate onto a regular z-grid as was on model levels before^M

    potemp(:,i_line-n_line(1)+1) = interp1 (z_slice(i_line-n_line(1)+1,:)/1e3,cont_slice(i_line-n_line(1)+1,:),zz(1).z); %interpolate onto a regular z-grid as was on model levels before
    u = interp1 (z_slice(i_line-n_line(1)+1,:)/1e3,u_slice(i_line-n_line(1)+1,:),zz(1).z);
    v = interp1 (z_slice(i_line-n_line(1)+1,:)/1e3,v_slice(i_line-n_line(1)+1,:),zz(1).z);
    %%%%%%%%%%%%%%%%%
    sp = sqrt( u.^2 + v.^2 );

    clear dir
    %                    for iuv=1:length(sp)
    theta2 = 180/pi * atan ( u ./ v );

    iuv=find(u==0 & v==0);
    %                        if u(iuv)==0 & v(iuv)==0
    dir(iuv) = 0;
    iuv=find(u>=0 & v>=0);
    %                        elseif u(iuv)>=0 & v(iuv)>=0
    dir(iuv) = theta2(iuv);
    iuv = find(u>0 & v<0);  %theta2 is negative
    %                        elseif u(iuv)>0 & v(iuv)<0  %theta2 is negative
    dir(iuv) = 180 + theta2(iuv);
    iuv = find(u<=0 & v<=0);
    %                        elseif u(iuv)<=0 & v(iuv)<=0
    dir(iuv) = 180 + theta2(iuv);
    iuv=find(u<0 & v>0);
    %                        elseif u(iuv)<0 & v(iuv)>0
    dir(iuv) = 360 + theta2(iuv); %theta2 is negative
    %                        end
    %                    end

    dir = dir.*pi/180; %convert to radians

    u_section(:,i_line-n_line(1)+1) = sp.*cos(pi/2-angle(i_line)-dir); %find component for horizontal axis of the cross section
    v_section(:,i_line-n_line(1)+1) = interp1 (z_slice(i_line-n_line(1)+1,:)/1e3,w_slice(i_line-n_line(1)+1,:),zz(1).z); %interpolate onto a regular z-grid as was on model levels before

end

%%%%%%%%%% make streamline %%%%%%%%%%%%
%%% terrain data is in terr_slice
hh=[0 4 11 21 33 45 51 55 57.5]/58.5; %relation between h and H0 for turning over streamlines taken from graph in Smith (1985)
HH=[3 4 5 6 7 8 8.5 8.75 9]*pi/6;

delL=[0 7 12.5 18 24 30.5 32.5 34.5 35]*3/87; %del*L values that correspond to the hh values


dx=1e3; %spacing in x
X3=[X(1):dx:X(end)];
del=zeros(size(X3));
terr3=interp1(X,terr_slice,X3);
terr_orig=terr3;
X_orig=X3;

ismooth=1;
if ismooth==1
    nfilter=20; bfilter=ones([1 nfilter])*1/nfilter;
    terr3=filter(bfilter,1,terr3);
    terr3(end-nfilter+1:end)=[]; X3(end-nfilter+1:end)=[];
    terr3(1:nfilter)=[]; X3(1:nfilter)=[];
    X_terr=X3;
end

%terr_options='effective mountain';
%set offset_eff=0 for no artificial reduction of mountain height 
        offset_eff = 0;
        terr2 = terr3 - offset_eff; %slice off the bottom of the mountain
        terr2(terr2<0)=0;

[terr_max imax]=max(terr2);        
        
H02=H0-offset_eff;
iH0=findheight_nearest(Y,H0);
H=Y(iH0);
N= sqrt( 9.81 ./ potemp(2:end) .* diff(potemp)./diff(Y) );
L= N(iH0-1)./u_section(iH0);

%%%%%%%%%%%% override L? %%%%%%%%%%%
L=0.0016;
L=0.002;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if H02*L<=max(HH)
    hmax=interp1(HH,hh,H02*L)/L; %find the maximum mountain allowed for given H0 according to relationship between h and H0 in Smith (1985)
else
    hmax=hh(end)/L; %approx 1/L
end

%hmax=100;
offset=terr_max - hmax;

%hmax=hmax*0.9;



if terr_max>hmax
    %
    iterate=7;
    switch iterate  %for case 0 leave offset as it is
        case 1   %iterate to get new H0 based changing hmax value
            offset_new=offset;
            go=1;
            while go==1
                offset=offset_new;
                H0_new=H02-offset_new;
                hmax_new=interp1(HH,hh,H0_new*L)/L
                offset_new=terr_max - hmax_new;
                if abs(offset_new-offset)<1; go=0; end
            end
            offset=offset_new;
            hmax=hmax_new;
            
            terr2=terr2-offset;
            terr2(terr2<0)=0; %produce effective mountain height
            [terr_max imax]=max(terr2); %new terr_max
            offset_terr=offset;
            
        case 2  %do one iteration to change the terrain but don't change H0                
            H0_new=H02-offset;
            hmax=interp1(HH,hh,H0_new*L)/L;
            offset_terr=terr_max - hmax;   
            
            terr2=terr2-offset_terr;
            terr2(terr2<0)=0; %produce effective mountain height
            [terr_max imax]=max(terr2); %new terr_max
            
        case 3       %change the terrain inline with the max possible mountain but don't change H0    
                     %i.e. always follow the max curve for a given H0
            offset_terr=terr_max - hmax;  
%            offset_terr=0;
            
            terr2=terr2-offset_terr;
            terr2(terr2<0)=0; %produce effective mountain height
            [terr_max imax]=max(terr2); %new terr_max   
            
            offset=0; %no offset of H0
            
        case 4       %change the terrain inline with the max possible mountain AND change H0    
                     %i.e. always follow the max curve for a given H0
            offset_terr=terr_max - hmax;  
%            offset_terr=0;
            
            terr2=terr2-offset_terr;
            terr2(terr2<0)=0; %produce effective mountain height
            [terr_max imax]=max(terr2); %new terr_max   
            
        case 5  %do one iteration to change the terrain and change H0                
            H0_new=H02-offset; %changes H0 to original H0-max mountain + 1/L
            hmax=interp1(HH,hh,H0_new*L)/L; %find new hmax for new H0

            offset=terr_max - hmax; %set so that newer still H0 is original H0-max mountain + new hmax 
%%%%%%%%%%%%%%%%            
%            offset = H02 - 6*pi/6/L + (2900-H0); %can override the value of H0 here (set to H02 - H0_required)
%%%            offset = 2900 - 8*pi/6/L; %can override the value of H0 here (set to H02 - H0_required)
            NL=9; offset = 1500-interp1(HH,hh,NL*pi/6)/L
%            hmax= 470.34;
%%%%%%%%%%%%%%%%            
            hmax=interp1(HH,hh,(H02-offset)*L)/L; %find new hmax for new H0
                        
         %   hmax=hmax*0.85;  
            
            offset_terr=terr_max - hmax;   %change terrain so that won't go above hmax
            
            terr2=terr2-offset_terr;
            terr2(terr2<0)=0; %produce effective mountain height
            [terr_max imax]=max(terr2); %new terr_max      
            
        case 6    %don't change the terrain inline with the max possible mountain but don't change H0    
                     %i.e. always follow the max curve for a given H0
            offset_terr=terr_max - hmax;  
%            offset_terr=0;
            
            terr2=terr2-offset_terr;
            terr2(terr2<0)=0; %produce effective mountain height
            [terr_max imax]=max(terr2); %new terr_max   
            
            offset=0; %no offset of H0
            
        case 7    %choose H0 as a function of different chosen hhat values  
            NL=4;    
            H0_choose=NL*pi/6;
            hmax=interp1(HH,hh,H0_choose)/L;
            offset = terr_max-hmax;
            
%            offset=0; %no offset of H0    
  
            H02=H0_choose/L+offset;
                        
            offset_terr=offset;  
%            offset_terr=0;
            
            terr2=terr2-offset_terr;
            terr2(terr2<0)=0; %produce effective mountain height
            [terr_max imax]=max(terr2); %new terr_max   
            
            
            
    end
    
    

    


else
    offset=0;
    offset_terr=0;
end


H03 = H02-offset; %H0 from the base of the max allowed mountain







i0=find(terr2>0);
inds=i0(end)+1:length(terr2);
terr2(inds)=terr2(inds)-(offset_terr-terr3(inds));

clear del
for i3=i0(1):length(terr2)
    [delA delB]=mountain_smith_solve_fun(terr2(i3),L,H03);
    
    if i3<imax  %if are on way up the mountain
        del(i3)=delA;
    else  %if coming back down again!
        if terr_max>hmax*0.95  %if mountain is big enough to reach apex of Smith's solution then
            del(i3)=delB;  %use the del value from the lower branch (transition to super-critical)
        else
            del(i3)=delA; %otherwise use the value for going back down the upper branch
        end
    end    
       
end

        inan=isnan(del);
        inan=find(inan==1);
        X3(inan)=[];
        del(inan)=[];
        
%         if length(inan)>=1
%             xnan=X3([inan(1)-1 inan(end)+1]);
%             delnan=del([inan(1)-1 inan(end)+1]);
%             del(inan)=interp1(xnan,delnan,X3(inan)); %linear interpolation (straight line between the last known values)
%         end
xstream=X3/1000;
ystream=(H0+del)/1000;

'Finished smith_hydraulic streamlines'









