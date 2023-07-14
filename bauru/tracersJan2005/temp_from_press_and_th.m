function [T,P]=temp_from_press_and_th(Grid,th,p,iz)
%function [T]=temp_from_press_and_th(Grid,th,p,iz)
%works for 3D,1D and single values of th and p at the moment (i.e. not 2D - easy to implement though)
%for single values need to give the index of the height in Grid
%note: th is the potential temperature perturbation

a=size(th);
idim=find(a>1); %the length of this gives the number of dimensions in th

if length(idim)==3 %3D
	RHOref = repmat(Grid.RHON,[1 a(1)-2 a(2)-2]);
	Pref = repmat(Grid.PREFN,[1 a(1)-2 a(2)-2]);         
	thref = repmat(Grid.THREF,[1 a(1)-2 a(2)-2]);
	
	th=permute(th(2:end-1,2:end-1,:),[3 1 2]);
    
elseif length(idim)==1 %1D vector
    RHOref = Grid.RHON;
    Pref = Grid.PREFN;
    thref = Grid.THREF;
    
elseif length(idim)==0   %single value
    RHOref = Grid.RHON(iz);
    Pref   = Grid.PREFN(iz);
    thref  = Grid.THREF(iz);
end
    
	
	P = squeeze(  p .* RHOref + Pref );
	
	T=(th+thref)./(1000e2./P).^0.286;