
function [table_vals] = ACSIS_Robson_paper_offline_SW_calcs_Sep2020_TABLE_full_vals(sw_trends,fdT_P1,fdT_P2,fscale)


%Full SW trend values from the model
    table_vals.FullModelPA = fdT_P1*fscale*sw_trends.trend_dat_box{1,1}.coeffs(2);
    table_vals.FullModelPB = fdT_P2*fscale*sw_trends.trend_dat_box{1,2}.coeffs(2);
    table_vals.FullModelPAun = fdT_P1*fscale*sw_trends.trend_dat_box{1,1}.uncer_max;
    table_vals.FullModelPBun = fdT_P2*fscale*sw_trends.trend_dat_box{1,2}.uncer_max;
    
    %Estimated SW trend values from the offline calc
    table_vals.FullCalcPA = fdT_P1*fscale*sw_trends.trend_dat_box_obs2{1,1}.coeffs(2);
    table_vals.FullCalcPB = fdT_P2*fscale*sw_trends.trend_dat_box_obs2{1,2}.coeffs(2);
    table_vals.FullCalcPAun = fdT_P1*fscale*sw_trends.trend_dat_box_obs2{1,1}.uncer_max;
    table_vals.FullCalcPBun = fdT_P2*fscale*sw_trends.trend_dat_box_obs2{1,2}.uncer_max;