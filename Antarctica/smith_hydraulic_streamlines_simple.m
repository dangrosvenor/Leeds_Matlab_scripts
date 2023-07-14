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
X3_orig=X3;
del=zeros(size(X3));

terr3=interp1(X,terr_slice,X3);
terr_orig=terr3;

terrain_type='AP';
%terrain_type='Witch';
terrain_type='heff';
%terrain_type='heff_2';

switch terrain_type
    case 'AP'
        
        X_orig=X3;

        ismooth=0;
        if ismooth==1
            nfilter=20; bfilter=ones([1 nfilter])*1/nfilter;
            terr3=filter(bfilter,1,terr3);
            terr3(end-nfilter+1:end)=[]; X3(end-nfilter+1:end)=[];
            terr3(1:nfilter)=[]; X3(1:nfilter)=[];
            X_terr=X3;
        end
        
    case 'Witch'
        a=10e3;
        h=1400;
        terr3=h*a^2./((X3-425e3).^2+a^2);
        
    case 'heff'
        terr3=heff*1000;
        
    case 'heff_2'
        NL=4.86;    
        H0_choose=NL*pi/6;  %H0_hat 
            
        hmax = mountain_smith_solve_hhat_for_given_critical_H0_hat(H0_choose)/L;
        imask=find(heff*1000<max(heff)*1000-hmax);        
%        offset = terr_max-hmax;
        
        terr3=heff*1000;    
        terr3(imask)=max(heff)*1000-hmax;
end

%terr_options='effective mountain';
%set offset_eff=0 for no artificial reduction of mountain height 
        offset_eff = 0;
        terr2 = terr3 - offset_eff; %slice off the bottom of the mountain
        terr2(terr2<0)=0;

[terr_max imax]=max(terr2);        
        
N= sqrt( 9.81 ./ potemp(2:end) .* diff(potemp)./diff(Y) );
L= N(iH0-1)./u_section(iH0);

%%%%%%%%%%%% override L? %%%%%%%%%%%
L=0.0016;
L=0.002;
%L=0.01/10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if terr_max>hmax
    %
    iterate=7;
    switch iterate  %for case 0 leave offset as it is
            
        case 7    %choose H0 as a function of different chosen hhat values  
%            NL=5.53;    
            H0_choose=NL*pi/6;  %H0_hat
            
%            hmax=interp1(HH,hh,H0_choose)/L;
            hmax = hhat_peak/L;    
            
            hmax = mountain_smith_solve_hhat_for_given_critical_H0_hat(H0_choose)/L;
            offset = terr_max-hmax;
            
%            offset=0; %no offset of H0    
  
            H02=H0_choose/L+offset;
                        
            offset_terr=offset;  
%            offset_terr=0;

%            [del_Houghton,H0_Houghton]=Houghton_find_thi_streamline_from_mountain_func_version(terr2,9,0.7);
            
            terr2=terr2-offset_terr;
            terr2(terr2<0)=0; %produce effective mountain height
            [terr_max imax]=max(terr2); %new terr_max   
            
             [del_Houghton,H0_Houghton]=Houghton_find_thi_streamline_from_mountain_func_version(terr2,9,0.7);           
            
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
        terr3(inan)=[];
        terr2(inan)=[];
        
%         if length(inan)>=1
%             xnan=X3([inan(1)-1 inan(end)+1]);
%             delnan=del([inan(1)-1 inan(end)+1]);
%             del(inan)=interp1(xnan,delnan,X3(inan)); %linear interpolation (straight line between the last known values)
%         end
xstream=X3/1000;
ystream=(H0+del)/1000;

thi=H03+del-terr2;
F= ( (1-cos(L*thi)) ./ (1 - L*del.*sin(L.*thi)) ).^(-2);


'Finished smith_hydraulic streamlines'









