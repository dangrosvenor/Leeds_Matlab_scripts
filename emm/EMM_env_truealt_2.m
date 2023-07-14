% do environmental profile for EMM from LEM data in GridDan(1)

emm_case='NewMexico_Stewart';

switch emm_case
case 'NewMexico_Stewart'
    %run soundings with New Meciso case selected (sound='newmexico')
    soundings
    
    fid=fopen('c:/documents and settings/login/my documents/emm_0012/TEXT/ENV_NewMexico_Stewart.txt','wt');
    
    nmax=length(heights); %max no. points for profiles
    maxh=heights(end)/1000; %km
        
    epsilon=0.622;
    
    dh=maxh/(nmax-1);
    
    %set up linearly space grid
    h(1)=heights(1);
    h(nmax)=maxh;
    h(2:nmax-1)=[1:nmax-2]*dh;
        
    fprintf(fid,'%d\n',nmax);
    for i=1:nmax
        p(i)=press(i);
        t(i)=temp(i);
        
        %now work out the dew point temperature from qv - temperature at which qv would be saturated
            e(i)=p(i)*qvap(i)/(qvap(i)+epsilon); %vapour pressure
        %now need to find temp for es corresponding to e(i)        
            Td(i)=fzero(@eminus,273,[],e(i)); %eminus just does SatVapPress - e(i) so that fzero can find the solution we require. 273 is just a place to start
        % fzero at
        
        % fprintf(fid,'%f %f %e %f\n',p(i)/100,t(i)-273.15,wg(i),h(i));
        fprintf(fid,'%f %f %f %f\n',press(i)/100,temp(i)-273.15,Td(i)-273.15,heights(i)/1000);
        
    end
    
case 'LEM_sounding'
    
    %fid=fopen('c:/cygwin/home/user/emm/ENV_input','wt');
    fid=fopen('c:/documents and settings/g/my documents/emm_0012/TEXT/ENV_input_5ppmv_1000km_3','wt');
    
    nmax=150; %max no. points for profiles
    maxh=21; %20,165m
    
    
    epsilon=0.622;
    
    dh=maxh/(nmax-1);
    
    %set up linearly space grid
    h(1)=0;
    h(nmax)=maxh;
    h(2:nmax-1)=[1:nmax-2]*dh;
    
    T=TempLES(GridDan(1)); %works out temperature from LEM potemp
    
    fprintf(fid,'%d\n',nmax);
    for i=1:nmax
        ii=findheight(GridDan(1).Z,h(i)*1e3); %index for height concerned
        p(i)=GridDan(1).PREFN(ii);
        t(i)=T(ii);
        
        %now work out the dew point temperature from qv - temperature at which qv would be saturated
        wg(i)=GridDan(1).OLQBAR(ii,1); %in kg/kg
        %    e(i)=p(i)*wg(i)/(wg(i)+epsilon); %vapour pressure
        %now need to find temp for es corresponding to e(i)
        
        %    Td(i)=fzero(@eminus,273,[],e(i)); %eminus just does SatVapPress - e(i) so that fzero can find the solution we require. 273 is just a place to start
        % fzero at
        
        %    fprintf(fid,'ps_Pa[%d] = %f; temp_K[%d] = %f; td_K[%d] = %f;\n',i-1,p(i)/100,i-1,t(i)-273.15,i-1,Td(i)-273.15);                
        fprintf(fid,'%f %f %e %f\n',p(i)/100,t(i)-273.15,wg(i),h(i));                
        
    end
    
end
                                   

fclose(fid);
'done emm environment'
    
    