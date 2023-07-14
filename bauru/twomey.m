

% Cccn=[1000:50:6000];
% Cccn=1000;
% Cccn=2000;
Cccn=2000;
Kccn=0.8;
Kccn=0.64
W=[0:.1:50]; 

%S=[0:0.1:2];

%W=50;

nd=Cccn.^(2/(Kccn+2.)) * ( 1.62e-3*(100*W).^1.5 /Kccn/beta(1.5,Kccn/2.) ).^(Kccn/(Kccn+2)); %Twomey paper and LES

%nd=Cccn^0.8*W.^0.3;
%nd=Cccn.*S.^Kccn;

%nd=Cccn*( 3.6*( 1.6e-3*(100*W).^1.5/Cccn ).^(1/(Kccn+2)) ).^Kccn; %EMM Twomey from peak supersat

S=(nd./Cccn).^(1./Kccn); %from N=Cccn*S^Kccn

%nd = 7.e-2 * (100.*W).^1.5;
%nd = 0.88e6 * (Cccn).^(2./(2.+ Kccn)) * nd.^(Kccn/(2.+ Kccn))  ; 

%nd=1e-6*0.88e6*Cccn.^(2/(Kccn+2))*(0.07*(100.*W).^1.5).^(Kccn/(Kccn+2)); %EMM Twomey number of primary droplets at cloud bsae



%figure;plot(S.*100,nd/1e6);
figure;
plot(W,S);
%plot(W,nd);


Cccn=500;
Kccn=0.8;
Kccn=0.64
W=[0:.1:50]; 

nd=Cccn.^(2/(Kccn+2.)) * ( 1.62e-3*(100*W).^1.5 /Kccn/beta(1.5,Kccn/2.) ).^(Kccn/(Kccn+2)); %Twomey paper and LES
S=(nd./Cccn).^(1./Kccn); %from N=Cccn*S^Kccn

hold on
plot(W,S,'r');



