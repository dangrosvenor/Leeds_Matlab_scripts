i3d=1;


%%% Notes on scripts to use to work out qdef values %%%%
% Run watervapourMay2005.m with graph=4876 and the correct height ranges selected
% Prob best to use the loadselect=[118 133] for 3D/ccn comparisons as these are the 150 km runs
% and so comparisons are fairer (was using 300/150 km comparisons before - prob best not change text
% aobut use of 300 km domain since have all the stuff about the 2d/3d comparisons
% Needed to chnage tables as the percentages in table 5 weren't consistent with those in dq table
% unless error could just be put down to accuracy of the quoted numbers? (i.e. used proper numbers in the
% calc of percentage increases but only put down to 2 or 3 sf in tables.) But there is still the problem
% that vapour above 16 km was used for percentages in 6.2 (based on Liu frequencies) rather than above 15.9 km
% But this is not related to dq or table 5 since they are for all heights bigger than trop and dq is for 15-16 or
% 16-17 km. So can just change this table on own

switch i3d
    %%%%%%%%%%%  3d
		case 1
            D=300e3;
            D=150e3;
		%L=1000e3; %length of dehyd zone (km)

	%	A = pi * D^2 / 4;   % area of dehydrated air produced by model
            A = D^2;
            qdef=4; %2d moistening
            qdef=0.9;
           % qdef=0.17;
            qdef=0.08;
            qdef=0.07;
            
%             qdef=0.03;
%             
             qdef=0.3;
             qdef=0.33; %high CCN vap 16-17 km
             qdef=1.16; %high CCN vap 15-16 km
             
%   %           qdef=0.092;
%            qdef=0.013;
% %             
%              qdef=1;
%  %            qdef=0.83;
  %           qdef=0.76;
  
%  qdef=0.7;

        qdef=0.2852; %3D case 16-17km vap 
%        qdef=0.3361; %3D CCN case 16-17km vap   17.9 percent
%        qdef=1.0791; %3D case 15-16km vap
%        qdef=1.1634; %3D CCN case 15-16km vap    7.81 percent   
        
% %%%%%%%%% aside - values for tot for table 5 for high CCN vs normal %%%%%%%%%%%%%   
%         qdef=0.3045; %3D case 16-17km tot 
%         qdef=0.4140; %3D CCN case 16-17km tot    36 percent
%         qdef=1.9431; %3D case 15-16km tot
%         qdef=4.0381; %3D CCN case 15-16km tot    108 percent   
        
            
		case 0    
         %%%%%%%%%%%  2d
            D=2000e3;
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2d 3rd dim factor here %%%%%%%%%%%%%%%%%%%%%%%%%  
		    D2=200e3;
		    A = D*D2;
		
			qdef=1;
			qdef=0.2;
%			qdef=0.1; %2d moistening


end

L=2000e3; %processing length
V = 5  * (3600*24*30);  % horiz processing speed
%dndt=24 / (1000e3)^2;  %rate of events per km2 / month

dndt = ((1.94+0.48+0.23+0.22+0.29))/100 * 10198  / (pi*(240e3)^2) / (51/30);  %no. 10 dbz echotop clouds reaching >16 km per area per month
%dndt = ((0.48+0.23+0.22+0.29))/100 * 10198  / (pi*(240e3)^2) / (51/30);  %no. 10 dbz echotop clouds reaching >17 km per area per month
dndt = ((0.23+0.22+0.29))/100 * 10198  / (pi*(240e3)^2) / (51/30);  %no. 10 dbz echotop clouds reaching >17 km per area per month

%dndt = (0.13+0.51)/100 * 782  / (pi*(240e3)^2) / (51/30);  %no. 40 dbz echotop clouds reaching 15-16 km per area per month
%dndt = (0.13)/100 * 782  / (pi*(240e3)^2) / (51/30);  %no. 40 dbz echotop clouds reaching 15-16 km per area per month
%dndt = (0.13 + 0.51+5.5+5.88+7.42)/100 * 782  / (pi*(240e3)^2) / (51/30);  %no. 40 dbz echotop clouds reaching 15-16 km per area per month

%dndt = (0.4)/100 * 2484  / (pi*(240e3)^2) / (51/30);  %no. 40 dbz echotop clouds reaching 15-16 km per area per month
%dndt = (0.4 + 0.52)/100 * 2484  / (pi*(240e3)^2) / (51/30);  %no. 40 dbz echotop clouds reaching 15-16 km per area per month

dndt=(0.29+0.22+0.23)/100*10198 / (pi*(240e3)^2) / (51/30); %10 dbz 3D + 3D high ccn 

%dndt=(0.52+0.4)/100*2484  / (pi*(240e3)^2) / (51/30);   %3D high CCN 35 dBZ
%dndt=(5.5+0.51+0.13)/100*782 / (pi*(240e3)^2) / (51/30);  %3D high CCN 40 dBZ

dQ = A * dndt * qdef * L / V  %calculation of total dq change for air parcel traversing dehydration zone
                              %(assuming no uplift) of length L at speed V km/month
                              
                              
                              
                              
% = 0.4 ppmv for qdef=1, L=1000, v=5 m/s, D=200 km, dndt = 24 /month over 1000x1000 km                              
%derived from:
% mean reduction of mixing ratio over an area A2 of MR qenv due to the presence of dehydrated air of 
%q ppmv covering area A is:
% dq = (A*q + (A2 - A)*qenv ) / A2  -  qenv
% dq = A*(q - qenv)/A2   define q - qenv = qdef

% total number of events per month over A2 area = dndt * A2 where dndt is rate of events occuring per unit area
% then rate over area A2 dqdt = A*qdef*dndt
% so total change in q of air experiencing this for time of T = L/V is
% dQ = A*qdef*dndt * L/V

%should really apply this interactively so that the proper qdef is calculated from evolving qenv
% *** no need to consider uplift !! *** but obviously only going to affect certain height range 
% so uplift will need to be considered after air has passed through the region
% e.g. if work out that 1 km layer reduced by 2 ppmv for reasonable passage time then
% uplift rate will determine how often air will need to experience such a processing through one of these
% regions before moist air starts to follow it into stratosphere.
