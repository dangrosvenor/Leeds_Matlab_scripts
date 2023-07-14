%mpace cloud top base and heights - and then also binned into regular lon
%bins

idat=0;
                
%                 dual=2;    
%                 xloc=[1 1];

                tol=0.1;                      

                ibase=find(mpace_height>-tol & mpace_height<tol);                
                cb_lons = mpace_lon;
                cb_lons(ihtot) = NaN;
                cb_lons=cb_lons(ibase);                                
                cloud_base_height = mpace_alt(ibase);
               

                itop=find(mpace_height>1-tol & mpace_height<1+tol);                
                ct_lons = mpace_lon;
                ct_lons(ihtot) = NaN;
                ct_lons=ct_lons(itop);
                cloud_top_height = mpace_alt(itop);
               