figname=['Microphysical Process Rates for two numrates codes'];
        
        xdat(1).x=sum(icediag_nums(1).i(:,:,col),2);
        %xdat(2).x=sum(icediag_nums(2).i(:,:,col),2);
        %xdat(3).x=sum(f*TotMassBudgetALL(GridDan(idir),sum(icediagALL(idir).i(:,t1:t2,[22 23 24]),3),GridDan(idir).t,iz,iz2),2)*300;
        %xdat(3).x=xdat(1).x(1:end)-xdat(2).x(1:end);
        
 z=GridDan(1).Z;       
        xlab='Process Rate (kg^{-1}s^{-1})';
        titlenam=figname;
        
             ydat(1).y=z/1000;
            labs(1).l=[dgs{col} 'for run 1'];
            
            ydat(2).y=z/1000;
            labs(2).l=[dgs{col} 'for run 2'];
            
            %ydat(3).y=z/1000;
            %labs(3).l='Advection';