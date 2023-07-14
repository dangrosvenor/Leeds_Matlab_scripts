%3d_radar - set iovride=1 and idontclose=0 in Diag_newt_3d
%

fnmin=1;
fnmax=44;
%fnmax=8;

idir=1;


for jjfile=fnmin:fnmax  %for1        
    fprintf(1,'jjfile=%d ',jjfile);
    if jjfile<10;
        fn=strcat('00',int2str(jjfile));
    elseif jjfile<100
        fn=strcat('0',int2str(jjfile));
    else
        fn=int2str(jjfile);
    end

	fname=['RUN0001.DG0' fn];
    
    FileName=[direcDan(idir).dir fname];
    
    ijustheader=1; %set this so just reads in header in call to diag_3d_newt 
    Diag_3d_newt; %want header so can get APARMS tosee if  have graupel volume
    ijustheader=0; %reset flag for later
    
    

ztot=0;
for it=1:4

	switch it
	case 4 %ice
        lab='Ice';
        im=6;
        in=7;
        iv=99;
	case 2 %snow
        lab='Snow';
        im=4;
        in=9;
        iv=99;
	case 3 %graupel
        lab='Graupel';
        im=5;
        in=8;
        if APARMS(IHS+66)==0
            iv=99;
        else
            iv=11; %if have graupel volume then is q=11 (tracer is q=10)
        end
            
	case 1 %rain
        lab='Rain';
        im=3;
        in=99;
        iv=99;
	end
    
    %   n=ThreeD.Q(ih1:ih2,:,:,2);
 %   q=ThreeD.Q(ih1:ih2,:,:,1);
 
%  if iovride==0
% 	in=1;
% 	nreads(in)=0; in=in+1; %U 1
% 	nreads(in)=0; in=in+1; %V 2
% 	nreads(in)=0; in=in+1; %W 3
% 	nreads(in)=0; in=in+1; %TH 4
% 	nreads(in)=0; in=in+1; %Q01 5
% 	nreads(in)=0; in=in+1; %Q02 6
% 	nreads(in)=0; in=in+1; %Q03 7
% 	nreads(in)=0; in=in+1; %Q04 8
% 	nreads(in)=0; in=in+1; %Q05 9
% 	nreads(in)=0; in=in+1; %Q06 10
% 	nreads(in)=0; in=in+1; %Q07 11
% 	nreads(in)=0; in=in+1; %Q08 12
% 	nreads(in)=0; in=in+1; %Q09 13
% 	nreads(in)=0; in=in+1; %Q10 14
% 	nreads(in)=0; in=in+1; %Q11 15
% 	nreads(in)=0; in=in+1; %Q12 16
% 	nreads(in)=0; in=in+1; %P (pressure pert/rho) 17
% end

clear nreads
nreads=zeros([1 17]);
nreads(4+im)=1;
nreads(4+in)=1;
nreads(4+iv)=1;

Diag_3d_newt; %load in new variables


    
 %    RHO=repmat(Grid.RHON,[1 size(ThreeD.Q,1) size(ThreeD.Q,2)]);  %using this approximation to rho (instead of 2d pressure and temp
 %   RHO=permute(RHO,[2 3 1]);                                                        % makes little difference (< 0.5 dbz when tested for typical field)
                                                             %Using this makes life easier for 3d runs (less mem needed)
            
%             [r,c,p]=size(TwoD.P);
%			 PRESS=permute(repmat(Grid.PREFN,[1 1 c]),[1 3 2])+TwoD.P;
		%	 Tempera=TwoD.TH2.*((TwoD.PP)./100000).^0.286;
			% 
			 % Calculate density
		%	 RHO=(TwoD.PP./287)./Tempera;
        izmin=1;
        izmax=length(Grid.Z);
        if it==1
            ztot=ztot+Radar_new_selectHM(Grid,ThreeD.Q(:,:,:,1),ThreeD.Q(:,:,:,1),izmin,izmax,it,[]);
        else
            if iv==99
                ztot=ztot+Radar_new_selectHM(Grid,ThreeD.Q(:,:,:,2),ThreeD.Q(:,:,:,1),izmin,izmax,it,[]);
            else
                ztot=ztot+Radar_new_selectHM(Grid,ThreeD.Q(:,:,:,2),ThreeD.Q(:,:,:,1),izmin,izmax,it,ThreeD.Q(:,:,:,3));
            end

        end            
          %  clear ThreeD RHO
          %  pack
            
end
    

rads=[30 20 15 10 40 35];  %was [30 20 15 10] until 8th Jan, 2007   added 35 dbz 28th March 2007        

            jj=jjfile;
            for irad=1:length(rads)
				[arad]=find(10*log10(ztot)>=rads(irad)); %all points with >=n dBz
				[arad2]=find(10*log10(smooth3(ztot,'box',[1 1 17]))>=rads(irad)); %all points with >=n dBz                
                [ir,jr,kr]=ind2sub(size(ThreeD.Q(:,:,:,1)),arad); %convert to array indices from linear ones
				ihs=unique(kr); %all the different height indices with points > 10 dBz
				
				n10dbz(1).n(1:length(Grid.Z),irad,jj)=0;
				n10dbz(1).nsmooth(1:length(Grid.Z),irad,jj)=0;
				for iradar=1:length(ihs)
                    n10dbz(1).n(ihs(iradar),irad,jj)=length(find(kr==ihs(iradar))); %number of points > n dBz at each height index with > 10 dBz
                    n10dbz(1).nsmooth(ihs(iradar),irad,jj)=length(find(kr==ihs(iradar))); %number of points > n dBz at each height index with > 10 dBz                    
				end
            end
            
            [maxOverall imax]=maxALL(10*log10(ztot));
            n10dbz(1).max(jj)=maxOverall;
            n10dbz(1).hmax(jj)=Grid.Z(imax(3))+620;
            
            
            jj=jjfile;
            for irad=1:length(rads)
				[arad]=find(10*log10(smooth3(ztot,'box',[3 3 17]))>=rads(irad)); %all points with >=n dBz                
                [ir,jr,kr]=ind2sub(size(ThreeD.Q(:,:,:,1)),arad); %convert to array indices from linear ones
				ihs=unique(kr); %all the different height indices with points > 10 dBz
				
				n10dbz2(1).n(1:length(Grid.Z),irad,jj)=0;
				for iradar=1:length(ihs)
                    n10dbz2(1).n(ihs(iradar),irad,jj)=length(find(kr==ihs(iradar))); %number of points > n dBz at each height index with > 10 dBz
				end
            end
            
            [maxOverall imax]=maxALL(10*log10(smooth3(ztot,'box',[3 3 17])));
            n10dbz2(1).max(jj)=maxOverall;
            n10dbz2(1).hmax(jj)=Grid.Z(imax(3))+620;
            
            jc=1;
            if jjfile==fnmin
                sortheaderDGS;
                dgstrDan(1).dg=dgstr;
                TimeAvDan(1).DGAV=TimeAv.DGAV;
            end
            
            icediags_5thSept_2005_32; %do diags too


end  %jj      


