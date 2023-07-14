function [dat2]=expand_swath_lin_extrap(dat,sampling_along,sampling_across,re)
%expand the 5km L2 swath data fields to cover the full extent of the swath
%as for the 1km fields by (ad-hoc) linear extrapolation
%function [dat2]=expand_swath_lin_extrap(dat,sampling_along,sampling_across,re)
%uses re to get the size of the full swath

%add in the normal dat
dat2(sampling_across(1):sampling_across(2),sampling_along(1):sampling_along(2)) = dat;

%size(dat) is similar to [1346 2026]
%Y is the across dimension (1346, no. rows in swath)
ddat_Y1=dat(2,:) - dat(1,:);
ddat_Y2=dat(end,:) - dat(end-1,:);

%X is the along dimension (2026, no. rows in swath)
ddat_X1=dat(:,2) - dat(:,1);
ddat_X2=dat(:,end) - dat(:,end-1);

%Left indices
L_inds = double([sampling_along(1)-1:-1:1]);
%Right indices
R_inds = double([1:size(re,2)-sampling_along(2)]);
%Bottom indices
B_inds = double([sampling_across(1)-1:-1:1]);
%Top indices
T_inds = double([1:size(re,1)-sampling_across(2)]);

%LHS strip
LHS_inds = repmat(L_inds,[size(dat,1) 1]);
dat_left = repmat(dat(:,1),[1 size(LHS_inds,2)]);
ddat_left = repmat(ddat_X1,[1 size(LHS_inds,2)]);
BL_lhs_arr = dat_left - LHS_inds.*ddat_left;
dat2(sampling_across(1):sampling_across(2),1:sampling_along(1)-1) = BL_lhs_arr;

%RHS strip
RHS_inds = repmat(R_inds,[size(dat,1) 1]);
dat_right = repmat(dat(:,end),[1 size(RHS_inds,2)]);
ddat_right = repmat(ddat_X2,[1 size(RHS_inds,2)]);
BR_rhs_arr = dat_right + RHS_inds.*ddat_right;
dat2(sampling_across(1):sampling_across(2),sampling_along(2)+1:size(re,2)) = BR_rhs_arr;

%Bottom strip
BS_inds = repmat(B_inds,[size(dat,2) 1]);
dat_bott = repmat(dat(1,:),[size(BS_inds,2) 1]);
ddat_bott = repmat(ddat_Y1,[size(BS_inds,2) 1]);
BS_arr = dat_bott - BS_inds'.*ddat_bott;
dat2(1:sampling_across(1)-1,sampling_along(1):sampling_along(2)) = BS_arr;

%Top strip
TS_inds = repmat(T_inds,[size(dat,2) 1]);
dat_top = repmat(dat(end,:),[size(TS_inds,2) 1]);
ddat_top = repmat(ddat_Y2,[size(TS_inds,2) 1]);
TS_arr = dat_top + TS_inds'.*ddat_top;
dat2(sampling_across(2)+1:size(re,1),sampling_along(1):sampling_along(2)) = TS_arr;


%  *** The corners *** 

%Bottom right corner
BR_corner_arr = repmat(R_inds,[sampling_across(1)-1 1]);

ddat_X2_new = dat2(1:sampling_across(1)-1,sampling_along(2)) - dat2(1:sampling_across(1)-1,sampling_along(2)-1);
dx_repmat = repmat(ddat_X2_new,[1 size(BR_corner_arr,2)]);
dat_repmat = repmat(dat2(1:sampling_across(1)-1,sampling_along(2)),[1 size(BR_corner_arr,2)]);

BR_corner_dat = dat_repmat + BR_corner_arr.*dx_repmat;
dat2(1:sampling_across(1)-1,sampling_along(2)+1:size(re,2)) = BR_corner_dat;


%Top right corner
TR_corner_arr = repmat(R_inds,[size(re,1)-sampling_across(2) 1]);

ddat_X2_new = dat2(sampling_across(2)+1:size(re,1),sampling_along(2)) - dat2(sampling_across(2)+1:size(re,1),sampling_along(2)-1);
dat_repmat = repmat(dat2(sampling_across(2)+1:size(re,1),sampling_along(2)),[1 size(TR_corner_arr,2)]);
dx_repmat = repmat(ddat_X2_new,[1 size(TR_corner_arr,2)]);

TR_corner_dat = dat_repmat + TR_corner_arr.*dx_repmat;
dat2(sampling_across(2)+1:size(re,1),sampling_along(2)+1:size(re,2)) = TR_corner_dat;


%Top left corner
TL_corner_arr = repmat(L_inds,[size(re,1)-sampling_across(2) 1]);

ddat_X1_new = dat2(sampling_across(2)+1:size(re,1),2) - dat2(sampling_across(2)+1:size(re,1),1);
dx_repmat = repmat(ddat_X1_new,[1 size(TL_corner_arr,2)]);
dat_repmat = repmat(dat2(sampling_across(2)+1:size(re,1),sampling_along(1)),[1 size(TL_corner_arr,2)]);

TL_corner_dat = dat_repmat - TL_corner_arr.*dx_repmat;
dat2(sampling_across(2)+1:size(re,1),1:sampling_along(1)-1) = TL_corner_dat;


%Bottom left corner
BL_corner_arr = repmat(L_inds,[sampling_across(1)-1 1]);

ddat_X1_new = dat2(1:sampling_across(1)-1,2) - dat2(1:sampling_across(1)-1,1);
dx_repmat = repmat(ddat_X1_new,[1 size(BL_corner_arr,2)]);
dat_repmat = repmat(dat2(1:sampling_across(1)-1,sampling_along(1)),[1 size(TL_corner_arr,2)]);

BL_corner_dat = dat_repmat - BL_corner_arr.*dx_repmat;
dat2(1:sampling_across(1)-1,1:sampling_along(1)-1) = BL_corner_dat;