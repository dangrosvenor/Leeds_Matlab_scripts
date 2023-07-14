function LWC_ad=Adiabatic(Grid)

ClBase=input('Enter the height of cloud base (m)');
%ClLevel=input('Enter the level of interest (m)');

% Get some parameters
ClBase_ind=find(Grid.Z<=ClBase);
%ClLevel_ind=find(Grid.Z<=ClLevel);
Temp_ClBase=Grid.THREF(ClBase_ind(end)).*(Grid.PREFN(ClBase_ind(end))./100000).^0.286;
%Temp_ClLevel=Grid.THREF(ClLevel_ind(end)).*(Grid.PREFN(ClLevel_ind(end))./100000).^0.286;
Temp_ClLevel=Grid.THREF.*(Grid.PREFN./100000).^0.286;
% Saturation mixing ratios at those heights
QBase=LEWW(Temp_ClBase,Grid.PREFN(ClBase_ind(end)));
%Qlevel=LEWW(Temp_ClLevel,Grid.PREFN(ClLevel_ind(end)));
Qlevel=LEWW(Temp_ClLevel,Grid.PREFN);

%L/Lad
LWC_ad=QBase-Qlevel;
%LWC=Grid.OLQBAR(ClLevel_ind(end),2);
%LWC_ad./LWC
plot(LWC_ad,Grid.Z)