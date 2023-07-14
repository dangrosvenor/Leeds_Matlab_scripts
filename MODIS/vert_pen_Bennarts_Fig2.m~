filedir = '/home/disk/eos1/d.grosvenor/modis_work/vert_pen_paper/Bennartz_data/Bennartz2017_Fig2_data/';

cvals = [10:10:70];

fid = fopen([filedir 'Corners.csv']);
corners = textscan(fid,'%f%f\n','Delimiter',',');
fclose(fid);
xcorners=corners{1}; ycorners=corners{2};

x=[]; y=[]; z=[];
for i=1:length(cvals)
    fid = fopen([filedir num2str(cvals(i)) '_contour.csv']);
    dat = textscan(fid,'%f%f\n','Delimiter',',');
    fclose(fid);
    xdat{i}=dat{1}; ydat{i}=dat{2};
    
    x = cat(1,x,xdat{i}(:));
    y = cat(1,y,ydat{i}(:)); %N.B. - y is actually the distance from the top of the page...   
    z = cat( 1,z,cvals(i)*ones(size(xdat{i}(:))) );
    
end

%griddata(x,y,z,300,100)


%% Work out the width and x0 and y0 from the corners for normalisation
x0 = min(xcorners);
x1 = max(xcorners);
y0 = min(ycorners);
y1 = max(ycorners);

dx = x1-x0;
dy = y1-y0;

%% Normalise the data
xn = (x-x0)/dx;
yn = (y-y0)/dy;

%% calculate the actual values from x and y
H0=150; H1=400; 
N0=20; N1=1000; %log10 axis

N = 10.^ (xn * (log10(N1) - log10(N0)) + log10(N0) ); %since xn = ( log10(N) - log10(20) ) / (log10(1000) - log10(20))
H = H1 + yn*(H0-H1); %swap H0 and H1 from what wouuld expect these since y is measured from the top of the page

%% interpolate the dH values for a grid of N and H values based on the
%% contour data
H1d = [H0:1:H1];
N1d = 10.^([log10(N0):0.01:log10(N1)]);
N1d = [N0:1:N1];
%N1d = [N0:1:300];

[N2d,H2d] = meshgrid(N1d,H1d);

dH = griddata(N,H,z,N2d,H2d);


%% Calculate dtau values
CTT = 278;
CTP = 850e2;
fad=0.8;
k=0.8;

[dtau,tauc,cw]=MODIS_dtau_func_N_H_dH(H2d(:),N2d(:)*1e6,dH(:),CTT,k,CTP,fad);
dtau = reshape(dtau,size(H2d));
tauc = reshape(tauc,size(H2d));




%% Calculate adiabatic model re
W2d = 0.5 * fad *cw * H2d.^2;
[reff,H_model,k_model,Q_model,cw_model,tau_model]=MODIS_re_func_W_and_Nd(W2d,N2d*1e6,CTT,CTP,fad);
re_ret = ( 1 - (dtau ./ tauc) ).^(1/5) .* reff;

%% Do the line plots
figure; pcolor(N2d,H2d,dH); shading flat; colorbar; set(gca,'xscale','log');
figure; pcolor(N2d,H2d,dtau); shading flat; colorbar; set(gca,'xscale','log');
figure; pcolor(N2d,H2d,tauc); shading flat; colorbar; set(gca,'xscale','log');

%% 2D histograms
%Run this script to get the values for the 4 clouds in Platnick
MODIS_vert_pen_representative_tau_calc

%2D histos
Y_driver = dtau; %the data
X_driver = tauc; %the data
Ybins_DRIVER = [1:0.1:5]; ichoose_Ybins=1;
Xbins_DRIVER = [5:0.2:15]; ichoose_Xbins=1;
vert_pen_Bennartz_Fig2_2D_histo_dtau_vs_tau
hold on
plot(tauc_Plat,tau_star21,'ko-','linewidth',lwidth,'markerfacecolor','k');

%re_model - re_ret vs tau
Y_driver = (reff - re_ret)*1e6; %the data
X_driver = tauc; %the data
Ybins_DRIVER = [0:0.1:10]; ichoose_Ybins=1;
Xbins_DRIVER = [5:0.2:75]; ichoose_Xbins=1;
vert_pen_Bennartz_Fig2_2D_histo_dtau_vs_tau


dtau2=dtau;
dtau2(N2d>300)=NaN;
Y_driver = dtau2; %the data
X_driver = tauc; %the data
vert_pen_Bennartz_Fig2_2D_histo_dtau_vs_tau
hold on
plot(tauc_Plat,tau_star21,'ko-','linewidth',lwidth,'markerfacecolor','k');



% figure; plot(tauc(:),dtau(:),'kx');
% figure; plot(tauc(:),dtau(:),'kx'); set(gca,'xlim',[5 15]);
% 
% dtau2=dtau;
% dtau2(N2d>300)=NaN;
% figure; plot(tauc(:),dtau2(:),'kx'); set(gca,'xlim',[5 15]);
