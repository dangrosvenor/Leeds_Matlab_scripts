function [delA,delB]=Smith_solve_dual_layer_solve_del(Ha,Hb,hhat,idisp)
%function [delA,delB]=Smith_solve_dual_layer_solve_del(Ha,Hb,hhat,idisp)

del_start=5;
dx=-0.005;

delA=hhat-Ha -1;
delB=delA-Hb+Ha -1; %set up so that will enter the while loop at the start

while hhat>delA+Ha | delA>delB+(Hb-Ha)

    sign_fun=sign(Smith_dual_layer_solve_del_fun(del_start,Ha,Hb,hhat));
    sign_fun_old=sign_fun;
    x=del_start;
    while sign_fun==sign_fun_old
        x_start=x;
        x=x+dx;
        sign_fun=sign(Smith_dual_layer_solve_del_fun(x,Ha,Hb,hhat));
    end

    range=[min(x,x_start) max(x,x_start)];

    delA=fzero(@Smith_dual_layer_solve_del_fun,range,[],Ha,Hb,hhat);
    delB = - sqrt( delA.^2 + ( (delA - hhat) / (Ha + delA - hhat) ).^2 ); %eqn (23)

    del_start=x;
    
    if abs(delB) > 10 | abs(delA) > 10
        if idisp==1
            fprintf(1,'NO SOLUTION FOUND');
        end
        delA=NaN;
        delB=NaN;
        break
    end

end
