%choose a range of Nd and a start and end range of LWP to simlulate LWP
%adjustment (introduce positive LWP vs Nd relationship).
% This is the "true" dataset.
N=[40:5:40*7]*1e6;
W0=80e-3; %Start with 80 g/m2 cloud
Wend=102e-3; %Roughly what Antti did in the paper
dW=(Wend-W0)/(length(N)-1); %Increase linearly - they did log-linear, though
W=W0:dW:Wend;

%Calculate reff and tau from the given LWP and Nd values assuming T=283K,
%P=850 hPa.
[reff,H,k,Q,cw,tau]=MODIS_re_func_W_and_Nd(W,N,283,850e2,1,0.8);

Nper_sample = 200;
Npts = length(W)*Nper_sample; %Will have 200 points with random errors for each re value. At the moment I just keep the 
%tau value the same as in the true dataset

%Now introduce errors into re based on a normal distribution of
%errors centred on zero with std dev of err_std percent
err_std = 15; %std dev of the distribution of the errors
pd=makedist('Normal','mu',0,'sigma',err_std); %creates a normal distribution object
x=random(pd,Npts,1); %Samples Npts random points from the PDF.

j=1;
clear reff_err
for i=1:length(reff)
    reff_err(i,:)=reff(i)*(1+x(j:j+Nper_sample-1)/100);
    j = j + Nper_sample;
end
%replicate the true tau values
tau_err = repmat(tau,[size(reff_err,2) 1]);
tau_err = tau_err';

%Calculate N and LWP values based on the reff values with the errors added
[N2,H2,W2,k2,Q2,cw2]=MODIS_N_H_func(tau_err(:),reff_err(:),'',NaN,283,1,0.8); %N in per cm3 LWP in kg/m2

%Bin the data in lnN space
%lnW=log(W2);
lnN=log(N2);
dlnN=max(lnN)-min(lnN);
Nbins=min(lnN):dlnN/20:max(lnN);
for i=1:length(Nbins)-1
    ii=find(lnN>Nbins(i) & lnN<Nbins(i+1));
    meW(i)=mean(W2(ii));
end

figure('color','w');
plot(exp(Nbins(1:end-1)),1e3*(meW));
xlabel('N_d (cm^{-3})');
ylabel('LWP (g m^{-2})');
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca,'ylim',[50 110]);

hold on
plot(N/1e6,W*1e3,'r');

clear leg_str
leg_str{1} = ['Reff % error, sigma = ' num2str(err_std)];
leg_str{2} = 'Imposed LWP adjustment (basis dataset with no errors)';
legend(leg_str);


