function [P,T,TH]=press3d(ix,iy,iz,Grid,ThreeD)
%function [press]=press3d(ix,iy,iz,Grid,ThreeD)

TH=Grid.THREF(iz)+ThreeD.TH1(ix,iy,iz);
P=ThreeD.P(ix,iy,iz)*Grid.RHON(iz) + Grid.PREFN(iz);
T= TH/( (1000e2/P).^0.286 );