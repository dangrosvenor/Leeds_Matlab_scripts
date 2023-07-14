function [h]=mountain_smith_solve_fun_h_for_given_del(del,L,H0)
%h=mountain height
%H0=upstream streamline height
%L=N/U
% if h<0
%     dmax=-100/L;
% else
%     dmax=-1.5/L;
% end
% 
% n=100;
% dD=-0.005e3;

range=[-8/L 1.5/L]; %range over which the function should change sign (and thus the solution lie)

h=fzero(@mountain_smith_func_h_for_given_del,range,[],del,L,H0);  %solve for a given del (displacement)

%     sign0=sign(mountain_smith_func(0,h,L,H0) );
%     for d0=[dD:dD:dmax]
%         sign1 = sign(mountain_smith_func(d0,h,L,H0));
%         if sign1~=sign0
%             break
%         end       
%     end
    
%     if sign1~=sign0        
%         del=fzero(@mountain_smith_func,[d0 d0-dD],[],h,L,H0);        
%     else
%          del=NaN;
%     end
%     
%     if ~isnan(del)
%         dmax=-3.5/L;
%         sign0=sign(mountain_smith_func(del+dD,h,L,H0) );
%         for d0=[del+2*dD:dD:dmax]
%             sign1 = sign(mountain_smith_func(d0,h,L,H0));
%             if sign1~=sign0
%                 break
%             end
%         end

%         if sign1~=sign0
%             del2=fzero(@mountain_smith_func,[d0 d0-dD],[],h,L,H0);
%         else
%             del2=NaN;
%         end
% 
%     else
%         
%             del2=NaN;  %no second del value as mountain is too small to go past apex of solution curve
%     end
%     
%     if h<0
%         del2=del;
%         del=NaN;
%     end


%'done mountain wave solve'







