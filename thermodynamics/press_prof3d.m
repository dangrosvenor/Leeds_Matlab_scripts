function [P,T,TH]=press_prof3d(ix,iy,Grid,ThreeD)
%gives a profile of pressure and temperature from the 3d fields at ix,iy
%function [P,T]=press_prof3d(ix,iy,Grid,ThreeD)

for iz=1:length(Grid.Z)
    [P(iz),T(iz),TH(iz)]=press3d(ix,iy,iz,Grid,ThreeD);
end
