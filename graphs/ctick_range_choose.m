switch modis_data_plot
    case 'Max Nd, no CF screening minus CF screening'

        %
%        cbar_vals = [0:20:300 400 500 750 1000];
        cbar_vals = [-1:1:5 10:5:25];
        cbar_vals_for_text = cbar_vals;
        x_cbar_vals=linspace(0,1,length(cbar_vals));
        
        
        if icolormap_cf_grey==1
            %experimental additional CF grayscale colorbar on the end
%            cbar_vals = [0:25:300 400 500 750 1000 1e12 2e20 4e20 6e20 8e20];
            cbar_vals = [0:20:300 400 500 1e12 2e20 4e20 6e20 8e20];
            x_cbar_vals=linspace(0,1,length(cbar_vals));

%            cbar_vals_for_text = [0:25:300 400 500 750 1000 0.0 0.2 0.4 0.6 0.8];
            cbar_vals_for_text = [0:20:300 400 500 0.0 0.2 0.4 0.6 0.8];            
            i_cbar_val_new_scale_start = find(cbar_vals==1e12);
        end
        
    case 'Number of Droplets Joint Histogram mean from selected days'
        %
%        cbar_vals = [0:20:300 400 500 750 1000];
        cbar_vals = [0:20:300 400 500];
        cbar_vals_for_text = cbar_vals;
        x_cbar_vals=linspace(0,1,length(cbar_vals));
        
        
        if icolormap_cf_grey==1
            %experimental additional CF grayscale colorbar on the end
%            cbar_vals = [0:25:300 400 500 750 1000 1e12 2e20 4e20 6e20 8e20];
            cbar_vals = [0:20:300 400 500 1e12 2e20 4e20 6e20 8e20];
            x_cbar_vals=linspace(0,1,length(cbar_vals));

%            cbar_vals_for_text = [0:25:300 400 500 750 1000 0.0 0.2 0.4 0.6 0.8];
            cbar_vals_for_text = [0:20:300 400 500 0.0 0.2 0.4 0.6 0.8];            
            i_cbar_val_new_scale_start = find(cbar_vals==1e12);
        end
        
    case 'RWP GCM (grid-box mean)'
        %
%        cbar_vals = [0:20:300 400 500 750 1000];
        cbar_vals = [0:20:300 400 500];
        cbar_vals_for_text = cbar_vals;
        x_cbar_vals=linspace(0,1,length(cbar_vals));
        
        
        if icolormap_cf_grey==1
            %experimental additional CF grayscale colorbar on the end
%            cbar_vals = [0:25:300 400 500 750 1000 1e12 2e20 4e20 6e20 8e20];
            cbar_vals = [0:20:300 400 500 1e12 2e20 4e20 6e20 8e20];
            x_cbar_vals=linspace(0,1,length(cbar_vals));

%            cbar_vals_for_text = [0:25:300 400 500 750 1000 0.0 0.2 0.4 0.6 0.8];
            cbar_vals_for_text = [0:20:300 400 500 0.0 0.2 0.4 0.6 0.8];            
            i_cbar_val_new_scale_start = find(cbar_vals==1e12);
        end
        
        
        
        
    case {'Cloud Depth Joint Histogram mean from selected days','Cloud Depth timeseries3 mean from selected days'}       
        %    
%        cbar_vals = [0:100:600 750 1000 1500];
        cbar_vals = [0:50:750];        
        cbar_vals_for_text = cbar_vals;
        x_cbar_vals=linspace(0,1,length(cbar_vals));
        
         if icolormap_cf_grey==1
            %experimental additional CF grayscale colorbar on the end
%            cbar_vals = [0:100:600 750 1000 1500 1e12 2e20 4e20 6e20 8e20];
            cbar_vals = [0:50:750 1e12 2e20 4e20 6e20 8e20];            
            x_cbar_vals=linspace(0,1,length(cbar_vals));

%            cbar_vals_for_text = [0:100:600 750 1000 1500 0.0 0.2 0.4 0.6 0.8];
            cbar_vals_for_text = [0:50:750 0.0 0.2 0.4 0.6 0.8];            
            i_cbar_val_new_scale_start = find(cbar_vals==1e12);

         end
         
case {'LWP Joint Histogram mean from selected days'}       
        %    
%        cbar_vals = [0:100:600 750 1000 1500];
        cbar_vals = [0:20:300];        
        cbar_vals_for_text = cbar_vals;
        x_cbar_vals=linspace(0,1,length(cbar_vals));
        
         if icolormap_cf_grey==1
            %experimental additional CF grayscale colorbar on the end
%            cbar_vals = [0:100:600 750 1000 1500 1e12 2e20 4e20 6e20 8e20];
            cbar_vals = [0:20:300 1e12 2e20 4e20 6e20 8e20];            
            x_cbar_vals=linspace(0,1,length(cbar_vals));

%            cbar_vals_for_text = [0:100:600 750 1000 1500 0.0 0.2 0.4 0.6 0.8];
            cbar_vals_for_text = [0:20:300 0.0 0.2 0.4 0.6 0.8];            
            i_cbar_val_new_scale_start = find(cbar_vals==1e12);

         end         

    case {'Cloud-Sat warm rain precip rate','Total surface precip rate GCM'}
        %    
%        cbar_vals = [0:100:600 750 1000 1500];
        cbar_vals = [0:0.005:0.03 0.035:0.01:0.06 0.07:0.04:0.21];        
%        cbar_vals = [0:0.001:0.01 0.05:0.04:0.21];        
        
        cbar_vals_for_text = cbar_vals;
        x_cbar_vals=linspace(0,1,length(cbar_vals));
        
         if icolormap_cf_grey==1
            %experimental additional CF grayscale colorbar on the end
%            cbar_vals = [0:100:600 750 1000 1500 1e12 2e20 4e20 6e20 8e20];
            cbar_vals = [0:20:300 1e12 2e20 4e20 6e20 8e20];            
            x_cbar_vals=linspace(0,1,length(cbar_vals));

%            cbar_vals_for_text = [0:100:600 750 1000 1500 0.0 0.2 0.4 0.6 0.8];
            cbar_vals_for_text = [0:20:300 0.0 0.2 0.4 0.6 0.8];            
            i_cbar_val_new_scale_start = find(cbar_vals==1e12);

         end         

         
    otherwise  %the set for Nd.

        cbar_vals = [0:25:225 400 500];
        cbar_vals_for_text = cbar_vals;
        x_cbar_vals=linspace(0,1,length(cbar_vals));


        if icolormap_cf_grey==1
            %experimental additional CF grayscale colorbar on the end
            %            cbar_vals = [0:25:300 400 500 750 1000 1e12 2e20 4e20 6e20 8e20];
            cbar_vals = [0:20:300 400 500 1e12 2e20 4e20 6e20 8e20];
            x_cbar_vals=linspace(0,1,length(cbar_vals));

            %            cbar_vals_for_text = [0:25:300 400 500 750 1000 0.0 0.2 0.4 0.6 0.8];
            cbar_vals_for_text = [0:20:300 400 500 0.0 0.2 0.4 0.6 0.8];
            i_cbar_val_new_scale_start = find(cbar_vals==1e12);
        end


end