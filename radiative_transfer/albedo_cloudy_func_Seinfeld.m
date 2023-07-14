function [Ac,tau] = albedo_cloudy_func_Seinfeld(W,Nm3,i_Liu)
%W in kg/m2, N in per m3

Q=2;
rhow=1000;
fad=1;
k=0.8; %not used if i_Liu==1


CTT = 278;
P=850*1e2; %Convert to Pa. Using constant pressure here. Could vary according to the pressure also? Although
%the pressure dependence is not very strong
Parr=ones(size(CTT))*P; %
[cw]=fad*adlwcgm2_just_robs(CTT,Parr) * 1e-3;  % kg/m4

%Make sure there are no negative numbers since it causes imaginary parts
W(W<0)=0;
Nm3(Nm3<0)=0; 

if exist('i_Liu') & i_Liu==1
    %x=1/k.^(1/3); y=0; %sets k=0.8 with no dependence on Nd and L for testing to see if get the same as usual approach.       
    y=0.14; %From Liu, ERL, 2008: beta = x * (Nd/L)^y relationship for dispersion (see Mulcahy, 2018)    
                       % = x*(1/M)^y
    x=0.07./(1000.^y); %In Liu values are quoted in cgs units. Nd/L (=1/M) has units grams^-1 in cgs = 1/M_grams. Sub in M_grams = 1000*M_kg :-
    %x=0.07;            % beta = x * (1/(1000*M_kg)).^y = (1/1000).^y * x*(1/M_kg)^y
                        % = x' * (1/M_kg)^y   with x' = x*(1/1000).^y 
    B = 3.*Q./(4*rhow) .* cw.^(2/3+y) ./ ( x.* (3/(4*pi*rhow)).^(1/3) .* Nm3.^(y-1/3) );
    H = (2.*W./cw).^(0.5);
    tau = B .* H.^(5/3+y) ./ (5/3+y);
else
            
    B = 2*sqrt(10)*rhow.^2*cw.^(0.5).*(5/9).^(5/2) ./ (k.*pi.*Q.^3);
    tau = (Nm3.*W.^(5/2)./B).^(1/3);
end

Ac = tau ./(tau+7.7); % (Eqn. 24.38 of Seinfeld and Pandis)