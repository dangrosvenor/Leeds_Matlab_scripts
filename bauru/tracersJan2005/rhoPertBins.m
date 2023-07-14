clear a2

tref=repmat(GridDan(idir).THREF,[1 length(GridDan(idir).Y1)]);
            pref=repmat(GridDan(idir).PREFN,[1 length(GridDan(idir).Y1)]); %ref p
            
            
            T=TwoDDan(idir).TH1+tref; %tot potemp
            Tav=repmat(mean(T,2),[1 length(GridDan(idir).Y1)]); %mean T at this point in time
            
            P=TwoDDan(idir).PP; %tot P
            Pav=repmat(mean(P,2),[1 length(GridDan(idir).Y1)]); %mean P at this time
            
            T=T./(1e5./P).^0.286; %tot temp
            tref=tref./(1e5./pref).^0.286; %tot temp
            
            Tav=Tav./(1e5./Pav).^0.286; %tot temp
            
            rho=P.*28.97e-3/8.3144./T;
            rhoref=pref.*28.97e-3/8.3144./tref;
            %rhoref=Pav.*28.97e-3/8.3144./Tav;
           
            
            rhomoist=rho.*(1+TwoDDan(idir).Q(:,:,1))./(1+1.608*TwoDDan(idir).Q(:,:,1));
            
            rhopert=rhomoist-rhoref;
            

binsRHO=[-6e-3:1e-4:10e-3];

binsRHO=[-10e-3:1e-4:10e-3];

minpps=[4.6 4.8 5];

ikm=132;

figure
cols={'r','b','g'};
for i=1:length(minpps)
	inon2=find(totw(ikm,:)*f<minpps(i) ); % & rhopert(ikm,:)>0);
	dat2=rhopert(ikm,inon2);
	a2(i,:)=binner(dat2,binsRHO);
	%plot(binsRHO(2:end),a);
    plot(binsRHO(2:end),a2(i,:),cols{i});
    hold on
end