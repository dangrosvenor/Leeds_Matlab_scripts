%3d_radar
ztot=0;
for it=1:4

	switch it
	case 4 %ice
        lab='Ice';
        im=6;
        in=7;
	case 2 %snow
        lab='Snow';
        im=4;
        in=9;
	case 3 %graupel
        lab='Graupel';
        im=5;
        in=8;
	case 1 %rain
        lab='Rain';
        im=3;
        in=3;        
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
nreads(4+im)=1;
nreads(4+in)=1;

Diag_3d_newt; %load in new variables


    
     RHO=repmat(GridDan(1).RHON,[1 size(ThreeD.Q,2) size(ThreeD.Q,3)]);  %using this approximation to rho (instead of 2d pressure and temp
                                                             % makes little difference (< 0.5 dbz when tested for typical field)
                                                             %Using this makes life easier for 3d runs (less mem needed)
            
%             [r,c,p]=size(TwoD.P);
%			 PRESS=permute(repmat(Grid.PREFN,[1 1 c]),[1 3 2])+TwoD.P;
		%	 Tempera=TwoD.TH2.*((TwoD.PP)./100000).^0.286;
			% 
			 % Calculate density
		%	 RHO=(TwoD.PP./287)./Tempera;

            ztot=ztot+Radar_new_select(Grid,ThreeD.Q(:,:,:,2),ThreeD.Q(:,:,:,1),izmin,izmax,RHO,it);
            
          %  clear ThreeD RHO
          %  pack
            
end
    

            