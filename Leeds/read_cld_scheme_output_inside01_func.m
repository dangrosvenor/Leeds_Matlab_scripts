function [iremove, n, alphaL, aL, bs, Qn, qcl_on_bs, qt, qsl, CF] = read_cld_scheme_output_inside01_func(read_file,format_str);


fid=fopen(read_file,'rt'); %is quicker to use fscanf that dlmread
dat = textscan(fid,format_str);
fclose(fid);




dat_col01 = 13;

%chop off the last line as can be incomplete
nrec = -1;
for i=1:length(dat)
    nrec = max(nrec,length(dat{i}) );
end



%also chop out ones that are not at the end of the n iteration in the cloud
%scheme
n = dat{dat_col01}(1:nrec-1);
dn=diff(n);
iremove = find(dn==1);

for i=1:length(dat)
    dat{i} = dat{i}(1:nrec-1);    
    dat{i}(iremove) = ''; %remove if difference is 1 (i.e. are incrementing n). Only want to keep when n drops from >1 to 1 (i.e. negative dn)    
end



icol=0;
n = dat{dat_col01 + icol}(:); icol=icol+1; 
alphaL = dat{dat_col01 + icol}(:); icol=icol+1; 
aL = dat{dat_col01 + icol}(:); icol=icol+1; 
bs = dat{dat_col01 + icol}(:); icol=icol+1; 
Qn = dat{dat_col01 + icol}(:); icol=icol+1; 
qcl_on_bs = dat{dat_col01 + icol}(:); icol=icol+1; 
qt = dat{dat_col01 + icol}(:); icol=icol+1; 
qsl = dat{dat_col01 + icol}(:); icol=icol+1; 



      

CF=NaN*ones(size(Qn));

iqn = find(Qn <= -1.0);
CF(iqn) = 0.0;

iqn = find(Qn > -1.0 & Qn <= 0.0);
CF(iqn) = 0.5 * (1. + Qn(iqn)) .* (1. + Qn(iqn));

iqn = find(Qn > 0.0 & Qn <= 1.0); 
CF(iqn) = 1. - 0.5 * (1. - Qn(iqn)) .* (1. - Qn(iqn));

iqn = find(Qn > 1.0);
CF(iqn) = 1.0;




'';