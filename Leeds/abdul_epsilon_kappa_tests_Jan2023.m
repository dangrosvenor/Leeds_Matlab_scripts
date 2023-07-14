qi=[0.3 0.3 0.4]; %mixing ratios of components a,b and c.
%qi=[0.3 0.3 0.0]; %mixing ratios of components a,b and c.
%qi=[0.3 0.0 0.0]; %mixing ratios of components a,b and c.
%qi=[0.3 0.0 0.3]; %mixing ratios of components a,b and c.
qtot = sum(qi);
Mw=1;
Mi=[1 1 1];
rhoi=[1 1 1];
rhow=1;
mui = [2 3 0]; %component 3 is assumed insoluble.
mui = [2 2 0]; %component 3 is assumed insoluble.

%qi=[0.6 0.4]; mui = [2 0]; rhoi=[1 1]; Mi=[1 1];

tot_top_01 = 0;
tot_bot_01 = 0;
clear eps ki_01 ki_02
for i=1:length(qi)
   eps(i) = qi(i)/qtot; %mass fraction as assumed in UKCA (only for soluble)
    
   ki_01(i) = rhoi(i)*Mw*mui(i)*eps(i)/(rhow*Mi(i));         
   
   if mui(i)>0
       tot_top_01 = tot_top_01 + ki_01(i) * qi(i)/rhoi(i);
       
       tot_bot_01 = tot_bot_01 + qi(i)/rhoi(i);
   end       
    
end

kmean_01 = tot_top_01/tot_bot_01


tot_top_02 = 0;
tot_bot_02 = 0;
for i=1:length(qi)
   %eps = qi(i)/qtot; %mass fraction as assumed in UKCA (only for soluble)      
   ki_02(i) = rhoi(i)*Mw*mui(i)/(rhow*Mi(i));       
   
   %if mui(i)>0
       
       tot_top_02 = tot_top_02 + ki_02(i) * qi(i)/rhoi(i);
       tot_bot_02 = tot_bot_02 + qi(i)/rhoi(i);
   %end       
    
end
kmean_02 = tot_top_02/tot_bot_02


Vi = qi./rhoi;
Vf = Vi/sum(Vi)





