iload=65;

%comp='lacieLap';
comp='lap';

switch comp
case 'lap'
    dirroot='c:/documents and settings/g/my documents/emm_0012/output/';
case 'uni'
    dirroot='c:/cygwin/home/user/emm/out201/';
end

if iload==1
    
	fid=fopen([dirroot 'ENV'],'rt');
	rh_env=fscanf(fid,'%f',[5 Inf]);
	fclose(fid);
	
	
	fid3=fopen([dirroot 'drops'],'rt');
	D=fscanf(fid3,'%f',[4 Inf]); %tim,flowing,z_l,droplet_N
	fclose(fid3);
               
	fid2=fopen([dirroot 'N_w'],'rt');
	%N=fscanf(fid2,'%f',[9 Inf]);
	N=fscanf(fid2,'%f',[7 Inf]);
    fclose(fid2);

end

epsilon=0.622;
f=1e6*28.97/18; %for conversion to ppmv from mixing ratio

% pemm=
% rhemm=rh_env(3,:);
% esemm=SatVapPress(rh_env(2,:),'roger','liq',0,0); %sat formula used in EMM is from Rogers and Yau
% eemm=rhemm.*esemm;
% qemm=epsilon*eemm./(pemm*100-eemm);
% 
% hemm=rh_env(5,:);

%takes a table from the EMM of N*Nv where N is the number of output data and Nv is the total number
%of vertical points and so includes several different times (and hence passes through the altitude
%the following sorts this into DD(i).d structures where i is flowing (e.g. =0,1 for updraught =2 for SCR
%=3 for downdraught. The .d arrays are sorted so that there is a seperate index for time and height

clear DD NN
n2=4; %no. data points
for j=0:3 %flowing
    [a b]=find(D(2,:)==j); %find all points where flowing ==j
    
    if length(a)>0 %only do if there is some data for this flowing value
        d=D(:,b); %get data for this flowing value
        dz=d(3,2:end)-d(3,1:end-1); 
        p=find(dz<0); %points where the altiude drops back down
        
        nt=length(p); %the number of such points determines the number of complete altitude passes
        if nt>=1 %if there are any complete passes
            DD(j+1).d=reshape(d(:,1:p(1)*nt),[n2 nt p(1)]); %reshape the complete passes - data,time(mins),height
        else
            p(1)=0; %so that calc below for partial pass works
        end
        
        DD(j+1).d(1:n2,nt+1,1:size(d,2)-p(1)*nt)=d(:,p(1)*nt+1:end); %add the partial last pass
    end
    
end


n2=7;
for j=0:3 %flowing
   % [a b]=find(N(2,:)==j); %find all points where flowing ==1
    [a b]=find(N(3,:)==j); %find all points where flowing ==1
    if length(a)>0
        d=N(:,b(1:end-1)); 
%        dz=d(8,2:end)-d(8,1:end-1);
        dz=d(1,2:end)-d(1,1:end-1);

        p=find(dz<0); %points where the altiude drops back down
        
        nt=length(p); %so nt is number of altitude passes - i.e. number of output times
        
        if nt>=1 %if there are any complete passes
            NN(j+1).n=reshape(d(:,1:p(1)*nt),[n2 p(1) nt]);
        else
            p(1)=0;
        end
        NN(j+1).n(1:n2,1:size(d,2)-p(1)*nt,nt+1)=d(:,p(1)*nt+1:end);
    end
    
    
%     last=1;
%     for i=1:length(p)
%         DD(j).d(i,1:4,1:p(i)-last+1)=d(:,last:p(i));
%         last=p(i);
%     end
    
end

'done'
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    