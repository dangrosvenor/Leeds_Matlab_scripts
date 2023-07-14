function [Nd] = nenes_func(var_in)
% Inputs into var_in structure (e..g var_in.N_aero, etc.) :-
% N_aero - total number conc of aerosl
% sig - width of size dist, deafult = 2
% Dg - median diameter of aerosol size dist
% w - updraft speed in m/s (N.B. - later converted to cm/s for Nenes)
% density - aerosol density. E.g. ammonium sulphate = 1777 kg/m3 (default)
% T - temperature in K. Default = 283 K
% P - air pressure in Pa

%Convert all of the variable names in the input structure to actual names
%for ease of use
name_struc='var_in'; %The name of the structure
names = eval(['fieldnames(' name_struc ');']);
for i=1:length(names)
    eval_str = [names{i} ' = ' name_struc '.' names{i} ';'];
    eval(eval_str);
end


if ~exist('sig')
    sig = 2; %Lognormal distribution width
end
if ~exist('T')
    T=283; %Temperature in K
end
if ~exist('density')
    density=1777;
end


w = w*100; %convert to cm/s as used by Nenes scheme


%Make an 800 bin aerosol size distribution (lognormal)
Nbins=800;

Dlow=Dg/10/sig;
Dhigh=Dg*10*sig;
Dspacing=( log(Dhigh)-log(Dlow) )/Nbins; %log spacing

logD=[log(Dlow):Dspacing:log(Dhigh)]; %equally spaced in log space

Dmid=(logD(2:end) + logD(1:end-1) )/2;
Dmid=exp(Dmid);

%Make a lognormal aerosol distribution
nd=lognormal(Dmid,N_aero,sig,Dg); %dN/dlogD

%Calculate the critical supersat of each dry aerosol size since this is
%what the Nenes scheme operates on. The critical supersat can then account for composition if
%needed.
% N.B - the formula below agrees approximately with the Scrit_Ghan2011 for
% ammonium sulphate
%[S,Dcrit]=Scrit(T,fliplr(exp(logD)/2)); %fliplr flips matrix so that low supersats will be first
%S in %

%Need to use a slightly different one for comparison to Ghan (2011) since it assumes a hygroscopicity of
% 0.7
[S]=Scrit_Ghan2011(T,fliplr(exp(logD)/2)); %fliplr flips matrix so that low supersats will be first


% The number in each bin - sum of this gives total N
% nd is dN/dlnD and Dspacing is dlnD
NN=nd*Dspacing;

switch act_routine
    case 'Nenes'
        for i=1:length(w)
            %Find the max supersat using fzero
            sm=fzerodan(@nenes_with_density,[min(S) max(S)],[],10e4,T,w(i),NN,S,density);
            %Plug back in to get the other values if needed
            [I(i),sp(i),smax(i)]=nenes_with_density(sm,10e4,T,w(i),NN,S,density);

            % Activate all the aerosol up to smax
            %    is=findheight(S,smax(i));
            is=findheight_nearest(S,smax(i));
            Nd(i)=sum(NN(1:is));
        end

    case 'Nenes consistent with Abdul'  %Special case of Nenes to use same constants (es, Dv, etc.) as Abdul UM scheme
        for i=1:length(w)
            %Find the max supersat using fzero
            sm=fzerodan(@nenes_consistent_with_abdul,[min(S) max(S)],[],P,T,w(i),NN,S,density);
            %Plug back in to get the other values if needed
            [I(i),sp(i),smax(i)]=nenes_consistent_with_abdul(sm,P,T,w(i),NN,S,density);

            % Activate all the aerosol up to smax
            %    is=findheight(S,smax(i));
            is=findheight_nearest(S,smax(i));
            Nd(i)=sum(NN(1:is));
        end

end


