%plots all the lines onto which to plot the sounding - can skip this stage by loading
%a figure of one drawn (and saved) earlier
figure
plot_tephi(213,303,213,100e2);

%the following requires the saturation vapour mixing ratio (in kg/kg), which can be calculated
%using the function q_sat:-

qsat=q_sat(tdat,pdat);

%plot the sounding
plot_tephi_data2(tdat,pdat,qdat,qsat,zdat,1); %the last flag can be 0,1 or 2. 
%when it is zero the script calcultes the CAPE, etc., but it might not work in all cicrumstances
%To switch this off set it to 1 or 2. These are the same except that 1 draws solid and 2 draws
%dotted lines.
%tdat - temperature sounding in K
%pdat - pressure in Pa
%qdat - vapour mixing ratio in kg/kg
%psat - saturation mixing ratio in kg/kg
%zdat - height in m