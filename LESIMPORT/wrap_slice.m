% wrap 2d  slice so that x,y=0 of domain is in the centre
%first dimesion of slice becomes the y-axis for the pcolor plot

clear slice slice2


idirs=[3];  %comment this if are running for all cases in Allimp...

iTwoD=1;

f=1e6*28.97/18;

iplot_wrap=0;
orient='horiz'; %horizontal slices
orient='vert';   %this only contains one case so unlikely to be required
orient='vert2';  %probably the one to use for vertical slices
%orient='min';

for idiri=1:length(idirs)
    idir=idirs(idiri);
   
    
switch orient
case 'horiz'
    
% ih=143; %134=16.5km 202=25km
ih=138;
 ih=135;
 ih=130;


% %ih=20;

%ih=142;
%ih=155;
%ih=146; %17.9 km
% ih=findheight((GridDan(1).Z+620)/1000 , 15.5);
% ih=findheight((GridDan(1).Z+620)/1000 , 14.8);
% ih=findheight((GridDan(1).Z+620)/1000 , 16.6);
 %ih=findheight((GridDan(1).Z+620)/1000 , 17.8);
 ih=findheight((GridDan(1).Z)/1000 , 6.8);
 ih=findheight((GridDan(1).Z)/1000 , 8);


(GridDan(idir).Z(ih)+620)/1000


field='w';
%field='vapour';
field='totice';
%field='v';
%field='tot_water';
%field='ice no';
%field='ice supersat';
%field='potemp';
%field='Tpert';
%field='ice supersat2';
%field='ice supersat3'; %using pressure_hslice etc.
%field='Q-field';


switch field
case 'Q-field'
	slice=squeeze(ThreeD.Q(2:end-1,2:end-1,ih,5));
    
case 'vapour'
	slice=f*squeeze(ThreeD.Q(2:end-1,2:end-1,ih));
    
case 'totcond'
    idir=2;
	slice=squeeze( sum(dat.Q(2:end-1,:,ih,1:4),4) );    
    
case 'potemp'

    
    meanprofile='yes';
    meanprofile='no';    
    switch meanprofile
    case 'yes'
        for ih=1:250
            
            
            th1=squeeze(ThreeD.TH1(2:end-1,2:end-1,ih));
            [ix iy]=size(th1);
%             ref = repmat(GridDan(idir).THREF(ih),[ix iy]);   
%             pot=th1+ref;
%             med=median(median(pot));
%             slice=pot-med;
%             me(ih)=mean(mean(slice));
            
            RHOref = repmat(GridDan(idir).RHON(ih),[ix iy]);
            Pref = repmat(GridDan(idir).PREFN(ih),[ix iy]);         
    
            P = squeeze(  ThreeD.P(2:end-1,2:end-1,ih) ) .* RHOref' + Pref';
          
              pot0 = repmat(GridDan(idir).THREF(ih),[ix iy])';
              T0 = pot0' ./ (1000e2./Pref).^0.286;
	      
              pot = pot0 + squeeze(  ThreeD.TH1(2:end-1,2:end-1,ih) );
              T = pot' ./ (1000e2./P').^0.286;
          
              me(ih) = mean(mean((T-median(median(T)) )));
            
            
        end
    end
    
	th1=squeeze(ThreeD.TH1(2:end-1,2:end-1,ih));
    [ix iy]=size(th1);
    ref = repmat(GridDan(idir).THREF(ih),[ix iy]);   
    pot=th1+ref;
    med=median(median(pot));
    slice=pot-med;
        
    
case 'Tpert'
	      %temp pert
          [ix iz]=size(squeeze(  ThreeD.TH1(2:end-1,2:end-1,ih) ));
          
          RHOref = repmat(GridDan(idir).RHON(ih),[iz ix]);
          Pref = repmat(GridDan(idir).PREFN(ih),[iz ix]);         
    
          P = squeeze(  ThreeD.P(2:end-1,2:end-1,ih) ) .* RHOref' + Pref';
          
          pot0 = repmat(GridDan(idir).THREF(ih),[iz ix])';
          
          Pmean=mean(mean(P));
          T0 = pot0' ./ (1000e2./Pmean).^0.286;
                    
%           P0=GridDan(idir).PREFN(ih);
%         
%         area=1;
         T0=icediagsALL(idir).i(ih,1,246)./area./(1e5./Pmean).^0.286;
        
	      
          pot = pot0 + squeeze(  ThreeD.TH1(2:end-1,2:end-1,ih) );
          T = pot' ./ (1000e2./P').^0.286;
          
          slice = (T-T0);
    
case 'w'
    slice=squeeze(ThreeD.W(2:end-1,2:end-1,ih));

case 'u'
    slice=squeeze(ThreeD.U(2:end-1,:,ih));
case 'v'
    slice=squeeze(ThreeD.V(2:end-1,:,ih));    
case 'tot_water'
    slice=f*squeeze(tot_water44(2:end-1,:,ih));   
case 'ice no'
    slice=squeeze(Q07(2:end-1,:,ih)); 
case 'ice supersat'
    slice=squeeze(TH1(2:end-1,:,ih));
    [ix iy]=size(slice);
    pot = repmat(GridDan(idir).THREF(ih),[ix iy]);
    slice=slice+pot; %overall theta
    pbar = repmat(GridDan(idir).PREFN(ih),[ix iy]);
    
    P = GridDan(idir).RHON(ih) * squeeze(ThreeD.P(2:end-1,:,ih)) + pbar  ;  %P is the pressure perturbation over the mean density
    T=slice./(1e5./P).^0.286;
    qsi=satvapPress(T,'lem','ice',P,1)/f; %satvappress gives in ppmv if 5th argument=1
    
    f=1e6*28.97/18;        
    vap=squeeze(Q01(2:end-1,:,ih));        
    slice=100*(vap-qsi)./qsi;
    
case 'ice supersat2'
    slice=squeeze(si(:,:,it));
    
case 'ice supersat3'
    it=44; %time to use slice at 17.9 km from 
    
    slice=potemp_hslice(idir).dat(2:end-1,2:end-1,it);
    [ix iy]=size(slice);
    pot = repmat(GridDan(idir).THREF(ih),[ix iy]);
    slice=slice+pot; %overall theta
    pbar = repmat(GridDan(idir).PREFN(ih),[ix iy]);
    
    P = GridDan(idir).RHON(ih) * squeeze(pressure_hslice(idir).dat(2:end-1,2:end-1,it)) + pbar  ;  %P is the pressure perturbation over the mean density
    T=slice./(1e5./P).^0.286;
    qsi=satvapPress(T,'lem','ice',P,1)/f; %satvappress gives in ppmv if 5th argument=1
    
    f=1e6*28.97/18;        
    vap=squeeze(vap_hslice(idir).dat(2:end-1,2:end-1,it));        
    slice=100*(vap-qsi)./qsi;    
    
    
%    slice=T;
    
end


[ix iy]=size(slice);

ixh=floor(ix/2);
iyh=floor(iy/2); %half way points


iwrap=0;   %N.B. there are two of these iwraps
if iwrap==1
    
    slice2(1:ix-ixh,:)=slice(ixh+1:ix,:);
    slice2(ix-ixh+1:ix,:) = slice(1:ixh,:);
    
    clear slice
    ixh=iyh;
    ix=iy;
    slice(:,1:ix-ixh)=slice2(:,ixh+1:ix);
    slice(:,ix-ixh+1:ix,:) = slice2(:,1:ixh);
    
end

qvmean=f*mean(mean(slice))


if iplot_wrap>=1
    figure
    if iplot_wrap==1
        pcolor(GridDan(idir).Y1(2:end-1)/1000,GridDan(idir).X1(2:end-1)/1000,slice);shading interp;colorbar  
    elseif iplot_wrap==2
        pcolor(slice);shading interp;colorbar  
    end
        
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'vert'
[ix iy iz]=size(ThreeD.Q);

iz=findheight(GridDan(1).Z+620,19.5e3);

izmin=120; %80
izs=[izmin:iz];

ixpos=ix-1;
ixpos=2;
ixpos=84;


ipos=findheight(GridDan(1).Y1/1000,2);
if ipos>=ixh+1;
    ixpos = ipos+1;
else
    ixpos = ix-ixh+ipos;
end



i=0;
%for ixpos=[148:151 2:10];
     i=i+1;
	slice2=squeeze(ThreeD.Q(ixpos,:,izs));



[ix iz]=size(slice2);

ixh=floor(ix/2); %half way points

clear slice
slice(1:ix-ixh,:)=slice2(ixh+1:ix,:);
slice(ix-ixh+1:ix,:) = slice2(1:ixh,:);




if iplot_wrap==1
    savedirname='c:/documents and settings/g/my documents/logbook/vapour_paper/3d_slices_dump8/';
    figure
    pcolor(GridDan(1).Y1(1:end)/1000,GridDan(1).Z(izs)/1000+0.62,f*slice');shading interp;
    caxis([0 10]);
    colorbar
    xstr=num2str(GridDan(1).Y1(ixpos)/1000);
    title( [xstr ' km']);
    exname=[savedirname int2str(i) '-' xstr 'km.emf'];
 	%print(gcf,'-dmeta',exname);
   %  close(gcf);
end

%end


case 'vert2'

iz=findheight(GridDan(1).Z+620,19.5e3);

izmin=120; %80
izs=[izmin:iz];

izs=1:length(GridDan(idir).Z);

i=0;

ixposits=[148:152 1:10];

ixposits=152; %original position (marked in cross section of draft)
ixposits=77; %150km case

%ixposits=20;

ixposits=2;
%ixposits=84;

 for ixpos=ixposits
     
%     for ixpos=1
     i=i+1;

var='vapour';    
var='v';
%var='w';
%var='tracer';
%var='tot';
var='pot';
%var='ice';
var='Tpert'
%var='ice2d';
var='totcond';
%var='Qfield';
%var='Qpert';


switch var
	case 'vapour'
         %vapour
		 slice2=squeeze(ThreeD.Q(2:end-1,ixpos,izs,2));
         varstr='Q';
	
	case 'ice'
         %all ice
         %idir=2;
%          if idir==3
%              idat=1;
%          else
%              idat=idir;
%          end
		 slice2=squeeze(  sum( ThreeDDan(idir).Q(2:end-1,ixpos,izs,4:6),4 )  );
         
       %  slice2=squeeze(  sum( ThreeDDan(idir).Q(ixpos,2:end-1,izs,4:6),4 )  );
         
          varstr='Q'; 
          
     case 'Qfield'
		  slice2=squeeze(  sum( ThreeDDan(idir).Q(2:end-1,ixpos,izs,1),4 )  );               
          varstr='Q';      
          
     case 'Qpert'
          qinds=[1];
          qref = repmat( sum( GridDan(idir).OLQBAR(:,qinds),2 ) ,[1 length(GridDan(idir).X1)-2])';
		  slice2=squeeze(  sum( ThreeDDan(idir).Q(2:end-1,ixpos,izs,qinds),4 )  ) - qref ;  
          
          varstr='Q';           
         
     case 'totcond'
         %all ice
         %idir=2;
		 slice2=squeeze(  sum( ThreeDDan(idir).Q(2:end-1,ixpos,izs,2:6),4 )  );    
         varstr='Q'; 
         
    case 'ice2d'
         %all ice
		 slice2=squeeze(  sum( TwoDDan(idir).Q(2:end-1,izs,4:6),3 )  ); 
         iTwoD=1;
         varstr='';   
         
	case 'pot'     
         %potemp
          [ix iz]=size(squeeze(  ThreeD.TH1(2:end-1,ixpos,izs) ));         
	      pot = repmat(GridDan(idir).THREF(izs),[1 ix]);
     	  slice2 = pot' + squeeze(  ThreeD.TH1(2:end-1,ixpos,izs) );
          
          iTwoD=1;
          varstr='TH2';
          
	case 'w'     
         %potemp
	%      pot = repmat(GridDan(idir).THREF(izs),[1 ix-2]);
%     	 slice2 =  squeeze(  ThreeD.W(2:end-1,ixpos,izs) ); 
         slice2=squeeze(  ThreeDDan(idir).W(2:end-1,ixpos,izs)  );    
         varstr='W'; 
	case 'v'     
     	 slice2 =  squeeze(  ThreeD.V(2:end-1,ixpos,izs) );           
%     	 slice2 =  squeeze(  ThreeD.V(20,:,izs) );   
         
         slice2=squeeze(  ThreeDDan(idir).V(2:end-1,ixpos,izs)  );    
         varstr='V'; 
         
     case 'tracer'
     	 slice2 =  squeeze(  ThreeD.Q(2:end-1,ixpos,izs,1) );   
         
     case 'tot'
     	 slice2 =  f*squeeze(  tot_water44(2:end-1,ixpos,izs,1) );    
     case 'Tpert'     
         %temp pert
          [ix iz]=size(squeeze(  ThreeDDan(idir).TH1(2:end-1,ixpos,izs) ));
          
          RHOref = repmat(GridDan(idir).RHON(izs),[1 ix]);
          Pref = repmat(GridDan(idir).PREFN(izs),[1 ix]);         
    
          P = squeeze(  ThreeDDan(idir).P(2:end-1,ixpos,izs) ) .* RHOref' + Pref';
          
          pot0 = repmat(GridDan(idir).THREF(izs),[1 ix])';
%          T0 = pot0' ./ (1000e2./Pref).^0.286;
          

	      
          pot = pot0 + squeeze(  ThreeDDan(idir).TH1(2:end-1,ixpos,izs) );
          T = pot' ./ (1000e2./P').^0.286;
          
          T0= repmat(median(T,2),[1 ix]);
          P0= repmat(median(P,1),[ix 1]);
          
%           filename=['c:/documents and settings/login/my documents/logbook/vapour_paper/pics/Mirvette/Tref'];
%     
%     fid=fopen(filename,'w');
%     fprintf(fid,'%s %s\n','Height','Temperature (K)');
%     for izt0=1:length(T0)
%         fprintf(fid,'%f %f\n',GridDan(idir).Z(izt0)+620,T0(izt0));
%     end
%     fclose(fid);
        
          
          slice2 = (T-T0)';
          slice2=(P-P0)/100;
     %   slice2 = (T)';
        varstr='T';
        varstr='P';

 end    


[ix iz]=size(slice2);

ixh=floor(ix/2); %half way points

iwrap=0;
if iwrap==1
clear slice
	slice(1:ix-ixh,:)=slice2(ixh+1:ix,:);
	slice(ix-ixh+1:ix,:) = slice2(1:ixh,:);
	slice=slice';
else
    slice=slice2';
end

if iTwoD==1
    comm=['TwoDDan(idir).' varstr '(:,2:size(slice,2)+1)=slice;'];
    eval(comm);
    comm=['TwoDDan(idir).' varstr '(:,size(slice,2)+2)=0;'];
    eval(comm);
end




if iplot_wrap==1
	savedirname='c:/documents and settings/g/my documents/logbook/vapour_paper/3d_slices2_dump8/';

    figure
    pcolor(GridDan(idir).Y1(2:end-1)/1000,GridDan(idir).Z(izs)/1000+0.62,slice);shading interp;
    %caxis([0 10]);
    colorbar
    xstr=num2str(GridDan(1).Y1(ixpos)/1000);
    title( [xstr ' km']);
    exname=[savedirname int2str(i) '-' xstr 'km.emf'];
 	%print(gcf,'-dmeta',exname);
   %  close(gcf);
end


end



case 'min'
minq=f*squeeze(min(min(ThreeD.Q(2:end-1,:,:)) ) );
figure
plot(minq,(GridDan(1).Z+620)/1000);
set(gca,'ylim',[15 19]);

    
    
    
end

end %for idir=1:length(idirs)
'done'
