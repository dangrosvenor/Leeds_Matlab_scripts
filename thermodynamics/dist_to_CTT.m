function [distCTT,T_CB] = dist_to_CTT(T,P,RH,CTT)

if isnan(T)==1 | isnan(P)==1 | isnan(RH)==1 | isnan(CTT)==1
    distCTT=NaN;
    T_CB = NaN;
    return
    
end

if CTT>=T
    distCTT = 0;
    T_CB = T;
else


    if RH<100
        [dist_to_CB,T_CB] = dist_to_LCL(T,P,RH); % in metres
    else
        dist_to_CB = 0;
        T_CB = T;
    end

    if CTT>=T_CB
        distCTT=NaN; %is an incompatability between RH and CTT
        T_CB=NaN;
    else

        [moist_ad_1,moist_ad_2] = moist_ad_lapse_rate(T_CB,P);
        %moist_ad_2 is the moist adiabatic temperature without the
        %C-C approximation. Also using temperature varation of L.
        %moist_ad_1 uses the C-C relationship (with T variation of L)


        distCTT = dist_to_CB + (T_CB - CTT)/moist_ad_2; %in metres

    end

end
