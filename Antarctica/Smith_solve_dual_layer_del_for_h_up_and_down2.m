function [delA,delB,hx,delA2,delB2,hx2,U,U2]=Smith_solve_dual_layer_del_for_h_up_and_down2(Ha,Hb,hpeak,U0,idisp)
%function [delA,delB]=Smith_solve_dual_layer_solve_del(Ha,Hb,hhat,idisp)


h0=hpeak; %mountain height
n=100;

dmax=5; %maximum  del_hat to solve up to

delA(1)=0;
delB(1)=0; %will be zero for h(x)=0 - won't try to solve as the f(del) curve only just touches y=0

%find all del values on the upslope
i=1;
d0_interval=0.005;
n_int=4;
del_start=d0_interval*n_int;
for h=h0/n:h0/n:h0-h0/n  %solves for each point up the mountain (from h0 back down to 0)
    i=i+1
    del_start=del_start-d0_interval*n_int; %try to start close to the previous solution
    %otherwise start to jump to different solution areas that seem to be valid but are not what we want
        
    %try the delA=0, delB=0 solution first
%     fSmith=Smith_dual_layer_solve_del_fun(0,Ha,Hb,h);
%     if abs(fSmith)<1e-20 
%             delA(i) = 0;
%             delB(i) = 0;
%             hx(i)=h; %store the mountain height used
%             continue
%     end

    

    
    delA_temp=h-Ha -1;
    delB_temp=delA_temp-Hb+Ha -1; %set up so that will enter the while loop at the start
    delA_current=1e9;
    

    
    while h>delA_temp+Ha | delA_temp>delB_temp+(Hb-Ha)
        near_zero=0;
        fSmith_store = 1e99;
        
        %if either of these conditions are met then the solution
        %is not valid so keep searching for the next one (think comes from eqns. 25 & 26 of Smith and Sun)
        sign0=sign(Smith_dual_layer_solve_del_fun(del_start,Ha,Hb,h));
%        for d0=[del_start+d0_interval:d0_interval:dmax]
        for d0=[del_start+d0_interval:d0_interval:dmax]
            %        sign1 = sign(mountain_smith_func(d0,h,L,H0));
            fSmith=Smith_dual_layer_solve_del_fun(d0,Ha,Hb,h);
            sign1 = sign(fSmith);
            if sign1~=sign0 
                break
            end
            %        sign0=sign1;
            
            tol=0.1;
            if abs(fSmith)<tol
                %if have a previously stored value then compare to that
                %otherwise store this first value
                if (near_zero==1 & abs(fSmith)<fSmith_store) | near_zero==0
                    fSmith_store=fSmith;
                    delA_store=d0;
                    h_min_store=h;               
                end

                near_zero=1;
                
            elseif near_zero==1 %if are no longer within the tol then break and use the stored value
                break
            end
                

        end
        

        
        if near_zero==1            
            delA_temp = delA_store;
            h=h_min_store;
            delB_temp = - sqrt( delA_temp.^2 + ( (delA_temp - h) / (Ha + delA_temp - h) ).^2 ); %eqn (23)delB_temp = - sqrt( delA(i).^2 + ( (delA(i) - hx(i)) / (Ha + delA(i) - hx(i)) ).^2 ); %eqn (23)
            %N.B. need to set delA_temp and delB_temp here as they are used in the while loop      
            
            del_start = delA_store;
                        
        else
            %        if d0==dmax-d0_interval
            if sign1==sign0 %if have come to the end of the loop without finding two points surrounding a solution
                disp('WARNING - solution not found');
                break
            end
        
            delA_temp=fzero(@Smith_dual_layer_solve_del_fun,[d0-d0_interval d0],[],Ha,Hb,h);
            delB_temp = - sqrt( delA_temp.^2 + ( (delA_temp - h) / (Ha + delA_temp - h) ).^2 ); %eqn (23)
    
        

    %        if h<delA_temp+Ha & delA_temp<delB_temp+(Hb-Ha) & abs(delA_temp)<abs(delA_current)
    %            delA_current=delA_temp; %store the currently favoured solution
    %            delB_current=delB_temp;
    %        end      
    
                del_start = d0;

        
        end
        


            delA(i)=delA_temp;    
            delB(i)=delB_temp;
            hx(i)=h; %store the mountain height used
        
        

    end
    
    
    
end

%now solve for going down the other side of the mountain
% - looking for solutions heading out to the more negative del values than at the peak

delA2=0;
delB2=0;
hx2=0;
U=0;
U2=0;

irun=1;
if irun==1


%the del for the h values on the downwind side must lie between here and del*L=-Inf
N=900;  %
i=0;
clear d2 hx2
d0_interval=-d0_interval; %reverse these to go downwards towards a very low delA value
dmax=-dmax;

%for h=h0-h0/n:-h0/n:-N*h0  %solves for each point up the mountain - from h0 back down to 0 and
                    % then below 0 (negative mountain height as can have an asymmetric mountain)
%for h=h0-4/n:-4/n:-4  
 
for h=h0:-4/n:-4
        i=i+1    
    
%    del_start=0; %actaully think delA=0 at h=0 so will be negative - N.B. this doesn't work
%because delA is not quite zero (should be, but some numerical inaccuracy)
if i>1
    del_start=delA2(i-1); %use the previous solution and work down from there
else
    del_start=delA(end); %otherwise use the solution for at the peak 
%    del_start=0;
end
    
    delA2(i)=h-Ha -1;
    delB2(i)=delA2(i)-Hb+Ha -1; %set up so that will enter the while loop at the start
    d0=[]; %reset this for the breakout condition later
    
    while (h>delA2(i)+Ha | delA2(i)>delB2(i)+(Hb-Ha) ) & del_start+d0_interval>dmax 
        %if either of these first 2 conditions are met then the solution
        %is not valid so keep searching for the next one (think comes from eqns. 25 & 26 of Smith and Sun)
        %the last condition is so that we stop solving when reach dmax
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
    
    
    if length(d0>0)    %otherwise have reach dmax
        delA2(i)=fzero(@Smith_dual_layer_solve_del_fun,[d0 d0-d0_interval],[],Ha,Hb,h);
        delB2(i) = - sqrt( delA2(i).^2 + ( (delA2(i) - h) / (Ha + delA2(i) - h) ).^2 ); %eqn (23)
        hx2(i)=h; %store the mountain height used
    else
        delA2(i)=NaN;
        delB2(i)=NaN;
        hx2(i)=h;
    end
end


%now find U as a function of h(x) for z=0
%represents the maximum wind speed and shows how it varies as the flow goes over the crest
%and down the other side of the mountain
dh = 0; %
Ahat = delA2.*cos(Ha+delA2); %(2.14) of Smith
Bhat = delA2.*sin(Ha+delA2); %(2.15)
z=hx2 + dh;

z=Ha+delA2; %the height of streamline a
ddel_dz=-Ahat.*sin(z)+Bhat.*cos(z); %found by differentiation of (2.9)
%U0=6.67;

A2 = (delA2 - hx2)./(Ha+delA2-hx2); %from first 2 of (22) in Smith & Sun

ddel_dz=A2; %differential of del wrt z

U2 = U0 * (1 - ddel_dz); %from (2.2) of Smith and just above (2) of Smith
% & Sun

if ~exist('hx')
    hx=0;
end

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



end