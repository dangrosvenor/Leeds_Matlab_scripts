function [F,h_h0]=mountain_Houghton_solve_cr_0_line(H0,gd)


F(1)=1;
h_h0(1)=0;

x_min=mountain_Houghton_solve_U0_for_h_H0(0.7*H0,H0,gd);

%solving:- function [cr]=mountain_Houghton_solve_cr(U0,Hc,H0,gd)

np=80;
i=1;
for Hc=H0/np:H0/np:H0-H0/np
%for Hc=3500*0.925
    i=i+1;

    x=max([x_min 1.005*mountain_Houghton_solve_U0_for_h_H0(Hc,H0,gd)]); %start at the U value for the critical condition. Mulitply by 1.01
    %to make sure that we are not below the critical line
    %solving function y=mountain_Houghton_fun_U0_for_h_H0(U0,Hc,H0,gd)
    
    sign_fun=sign(mountain_Houghton_solve_cr(x,Hc,H0,gd));
    sign_fun_old=sign_fun;

    while sign_fun==sign_fun_old
        x_start=x;
        x=x+1;
        sign_fun=sign(mountain_Houghton_solve_cr(x,Hc,H0,gd));  %find the approx U values where cr becomes zero
    end

    range=[x_start x]; %U0 values

    U=fzero(@mountain_Houghton_solve_cr,range,[],Hc,H0,gd); %solve more precisely using the roughly calculated range above

    F(i) = U/sqrt(gd*H0);
    h_h0(i)=Hc/H0;
%    Hc/H0

end
