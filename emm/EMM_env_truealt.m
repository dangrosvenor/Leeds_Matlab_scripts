% do environmental profile for EMM from LEM data in GridDan(1)


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
                                   

fclose(fid);
'done emm environment'
    
    