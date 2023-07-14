function [delA,delB,hx,delA2,delB2,hx2,U,U2]=Smith_solve_dual_layer_del_for_h_up_and_down(Ha,Hb,hpeak,U0,idisp)
%function [delA,delB]=Smith_solve_dual_layer_solve_del(Ha,Hb,hhat,idisp)


h0=hpeak; %mountain height
n=100;

dmax=-5; %minimum  del_hat to solve down to

%find all del values on the upslope
i=0;
d0_interval=-0.02;
for h=h0:-h0/n:0  %solves for each point up the mountain (from h0 back down to 0)
    i=i+1;
    
    del_start=5; %delA can be positive dual layer solution
    
    delA(i)=h-Ha -1;
    delB(i)=delA(i)-Hb+Ha -1; %set up so that will enter the while loop at the start
    
    while h>delA(i)+Ha | delA(i)>delB(i)+(Hb-Ha) %if either of these conditions are met then the solution
        %is not valid so keep searching for the next one (think comes from eqns. 25 & 26 of Smith and Sun)
        sign0=sign(Smith_dual_layer_solve_del_fun(del_start,Ha,Hb,h));
        for d0=[del_start+d0_interval:d0_interval:dmax]
            %        sign1 = sign(mountain_smith_func(d0,h,L,H0));
            fSmith=Smith_dual_layer_solve_del_fun(d0,Ha,Hb,h);
            sign1 = sign(fSmith);
            if sign1~=sign0 | abs(fSmith)<1e-20
                break
            end
            %        sign0=sign1;

        end

        delA(i)=fzero(@Smith_dual_layer_solve_del_fun,[d0 d0-d0_interval],[],Ha,Hb,h);
        delB(i) = - sqrt( delA(i).^2 + ( (delA(i) - h) / (Ha + delA(i) - h) ).^2 ); %eqn (23)
        hx(i)=h; %store the mountain height used
        
        del_start = d0;
        
        

    end
    
end

%now solve for going down the other side of the mountain
% - looking for solutions heading out to the more negative del values than at the peak

%the del for the h values on the downwind side must lie between here and del*L=-Inf
N=3;  %
i=0;
clear d2 hx2

for h=h0-h0/n:-h0/n:-N*h0  %solves for each point up the mountain - from h0 back down to 0 and
                    % then below 0 (negative mountain height as can have an asymmetric mountain)
    i=i+1;    
    
%    del_start=0; %actaully think delA=0 at h=0 so will be negative - N.B. this doesn't work
%because delA is not quite zero (should be, but some numerical inaccuracy)
if i>1
    del_start=delA2(i-1); %use the previous solution and work down from there
else
    del_start=delA(1); %otherwise use the solution for at the peak 
end
    
    delA2(i)=h-Ha -1;
    delB2(i)=delA2(i)-Hb+Ha -1; %set up so that will enter the while loop at the start

    
    while h>delA2(i)+Ha | delA2(i)>delB2(i)+(Hb-Ha) %if either of these conditions are met then the solution
        %is not valid so keep searching for the next one (think comes from eqns. 25 & 26 of Smith and Sun)
        sign0=sign(Smith_dual_layer_solve_del_fun(del_start,Ha,Hb,h));
        for d0=[del_start+d0_interval:d0_interval:dmax]  %
            %        sign1 = sign(mountain_smith_func(d0,h,L,H0));
            sign1 = sign(Smith_dual_layer_solve_del_fun(d0,Ha,Hb,h));
            if sign1~=sign0
                break
            end
            %        sign0=sign1;

        end

        delA2(i)=fzero(@Smith_dual_layer_solve_del_fun,[d0 d0-d0_interval],[],Ha,Hb,h);
        delB2(i) = - sqrt( delA2(i).^2 + ( (delA2(i) - h) / (Ha + delA2(i) - h) ).^2 ); %eqn (23)
        hx2(i)=h; %store the mountain height used
        
        del_start = d0;

    end
    
    
    
    
    delA2(i)=fzero(@Smith_dual_layer_solve_del_fun,[d0 d0-d0_interval],[],Ha,Hb,h);
    delB2(i) = - sqrt( delA2(i).^2 + ( (delA2(i) - h) / (Ha + delA2(i) - h) ).^2 ); %eqn (23)
    hx2(i)=h; %store the mountain height used
end


%now find U as a function of h(x) for z=0
%represents the maximum wind speed and shows how it varies as the flow goes over the crest
%and down the other side of the mountain
dh = 0; %to find speeds along a line that is 0.1/L above the mountain surface
Ahat = delA2.*cos(Ha+delA2); %(2.14) of Smith
Bhat = delA2.*sin(Ha+delA2); %(2.15)
z=hx2 + dh;
z=Ha+delA2; %the height of streamline a
ddel_dz=-Ahat.*sin(z)+Bhat.*cos(z); %found by differentiation of (2.9)
%U0=6.67;

A2 = (delA2 - hx2)./(Ha+delA2-hx2);

ddel_dz=A2;

%U0=6.67;
U2 = U0 * (1 - ddel_dz); %from (2.2) of Smith


A = (delA - hx)./(Ha+delA-hx);

ddel_dz=A;

%U0=6.67;
U = U0 * (1 - ddel_dz); %from (2.2) of Smith





'done mountain wave solve'








% %%%%%
% del_start=5;
% dx=-0.005;
% 
% delA=hhat-Ha -1;
% delB=delA-Hb+Ha -1; %set up so that will enter the while loop at the start
% 
% while hhat>delA+Ha | delA>delB+(Hb-Ha)
% 
%     sign_fun=sign(Smith_dual_layer_solve_del_fun(del_start,Ha,Hb,hhat));
%     sign_fun_old=sign_fun;
%     x=del_start;
%     while sign_fun==sign_fun_old
%         x_start=x;
%         x=x+dx;
%         sign_fun=sign(Smith_dual_layer_solve_del_fun(x,Ha,Hb,hhat));
%     end
% 
%     range=[min(x,x_start) max(x,x_start)];
% 
%     delA=fzero(@Smith_dual_layer_solve_del_fun,range,[],Ha,Hb,hhat);
%     delB = - sqrt( delA.^2 + ( (delA - hhat) / (Ha + delA - hhat) ).^2 ); %eqn (23)
% 
%     del_start=x;
%     
%     if abs(delB) > 10 | abs(delA) > 10
%         if idisp==1
%             fprintf(1,'NO SOLUTION FOUND');
%         end
%         delA=NaN;
%         delB=NaN;
%         break
%     end
% 
% end
