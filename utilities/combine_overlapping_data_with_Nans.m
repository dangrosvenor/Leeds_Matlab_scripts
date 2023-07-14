function [dat_out]=combine_overlapping_data_with_Nans(dat_cell,i_combine)

N=0;
tot=0;
for is=1:length(i_combine)
    dat01 = dat_cell{i_combine(is)};
    %   is=19; dat02 = SW_TOA_up_CERES_reg{is};
    %i=find(isnan(dat01)==1 & isnan(dat02)==1);
    i01=find(isnan(dat01)==1);
    %   i02=find(isnan(dat02)==1);

    dat01(i01)=0;
    %   dat02(i02)=0;

    N01=ones(size(dat01)); N01(i01)=0;
    %   N02=ones(size(dat02)); N02(i02)=0;
    N=N+N01; %total number of non-NaN points
    tot=tot+dat01;


end

N(N==0)=NaN;
dat_out = tot./N;