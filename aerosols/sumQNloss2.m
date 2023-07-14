%run multipop,  Combined_Sdist4 & massbin

%multipop properties:-
	% T=276; %283
	% Nbins=1e4;
	% 
	% 
	% %average background nulcei mode
	% Dg=0.016e-6; %mode diameter
	% sig=1.7;
	% Ntot=6400e6; %total no. conc m^-3

%Combined_Sdist4 properties:-
	% addstartbins=0;
	%     
	% nsbins=33; %33 for 30 bins
	% logflag=0;
	% itmax=5000;
	% 
	% weight=1;
	% 
	% Sstart=min([Sc(1).s Sc(2).s]);
	% Send=max([Sc(1).s(end) Sc(2).s(end)]);
	% 
	% 
	% Send=5.5/100;

%massbin properties:-
	% tol=6.5e-17;
	% 
	% tol=2e-17; %2.5e-17
	% weight=0;
	% massbinflag=1;



clear summ sumn

dt(2:length(times))=times(2:end)-times(1:end-1);
dt(1)=dt(2);



dgs='QNloss';
dgfind=findhead(dgs,dgstrDan(1).dg);

mass=repmat(mmav2,[length(Grid.Z) 1]);

dz(2:end)=Grid.Z(2:end)-Grid.Z(1:end-1);
dz(1)=0;
dz1=repmat(dz,[length(mmav2) 1])';




for t=1:size(diag(1).dg,3)


    summ(:,:,t)=diag(1).dg(:,dgfind(2:end),t).*mass.*dt(t).*dz1;
    sumn(:,:,t)=diag(1).dg(:,dgfind(2:end),t).*dt(t).*dz1;

end

