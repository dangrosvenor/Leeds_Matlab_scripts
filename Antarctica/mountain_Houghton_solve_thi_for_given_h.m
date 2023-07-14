function thi=mountain_Houghton_solve_thi_for_given_h(hx,uA,hA,gd,side)
%function thi=mountain_Houghton_solve_thi_for_given_h(hx,uA,hA,gd,side)
%finds the fluid depth, thi, for the given mountain height at location x
%side is a flag to say which side of the mountain is required
%side='windward' for the windward side and 'lee' for the lee side.


    K4=uA*hA;
    thi_crit=K4^(2/3)*gd^(-1/3); %critical value at mountain crest
    
    %if know del at the mountain top then we know thi_crit and so u_crit from u_crit=sqrt(g*thi_crit)
    %then can calculate K4=u_crit*thi_crit, which gives hA, if we know uA, from uA*hA=K4
    %Want to try and get heff from the displacement of the streamline and also consider
    %the displacement at the mountain crest to calculated hA?

%choose a range for thi in which a solution lies
%range=[100 1500];
%range=[1500 15000];

switch side   %are two solutions, one for the lee and one for the windward side
    case 'windward'
        range=[thi_crit hA];
    case 'lee'
        range=[1 thi_crit];
end

    
[H0_crit]=mountain_Houghton_solve_H0_for_h(hx,uA,gd);

%if this is the mountain peak then H0_crit will match hA for the given supplied hx (the peak)
if abs(hA-H0_crit)<1;
    thi=thi_crit; %critical value at mountain crest
else
    thi=fzero(@mountain_Houghton_fun_thi_for_given_h,range,[],hx,uA,hA,gd);
end

