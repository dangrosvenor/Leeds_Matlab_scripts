clear lnb2d lnb2d_vap lnb2d_tot

add_ground_height=0.62; %height to add to the vertical axis to account for level of ground abv msl.
f=1e6*28.97/18;


totw=sum(TwoD.Q(:,:,[1:6]),3); %jc is the number of the eg. GridDan(jc) where the data is stored
vap=sum(TwoD.Q(:,:,[1]),3);   %j just loops from 1:n
            
%minpps=[3.67 5 1 2 3 4]/f;
minpps=[5]/f;

pref=repmat(GridDan(idir).PREFN,[1 length(GridDan(idir).Y1)]); %ref p

thref=repmat(GridDan(idir).THREF,[1 length(GridDan(idir).Y1)]); %potemp
tref=thref./(1e5./pref).^0.286; %convert to temp


                        
T=TwoD(idir).TH1+thref; %tot potemp
Tav=repmat(mean(T,2),[1 length(GridDan(idir).Y1)]); %mean pot temp at this point in time

P=TwoD(idir).PP; %tot P
Pav=repmat(mean(P,2),[1 length(GridDan(idir).Y1)]); %mean P at this time

Tav=Tav./(1e5./Pav).^0.286; %tot temp

%izmin=130;
%izmax=130;

for ipps=1:length(minpps)

    for ikm=izmin:izmax  %size(Tav,1)   %size(totw,1)
        ikm
        inon=find(vap(ikm,:)<minpps(ipps));
        inon2=find(totw(ikm,:)<minpps(ipps));
       % inon=1:size(Tav,2);
%         for iy=1:length(inon)
%                 [ttemp,thtemp,ptemp,lnb2d(ikm,iy)]=lnb(T(ikm,inon(iy)),P(ikm,inon(iy)),Pav(:,1),Tav(:,1),GridDan(1).Z);
%         end
		lnb2d(ikm,1:size(Tav,2))=NaN;
        lnb2d_vap(ikm,1:size(Tav,2))=NaN;
        lnb2d_tot(ikm,1:size(Tav,2))=NaN;
        if length(inon)>0
      %      [ttemp,thtemp,ptemp,lnb2d(ikm,1:length(inon))]=lnb(T(ikm,inon),P(ikm,inon),pref(:,1),tref(:,1),GridDan(1).Z);
       %     [ttemp,thtemp,ptemp,lnb2d(ikm,1:length(inon))]=lnb(T(ikm,inon),P(ikm,inon),Pav(:,1),Tav(:,1),GridDan(1).Z);
       
       %using av profs 
           %[ttemp,thtemp,ptemp,lnb2d(ikm,1:length(inon))]=lnb_ice2(T(ikm,inon),P(ikm,inon),Pav(:,1),Tav(:,1),GridDan(1).Z,squeeze(TwoD.Q(ikm,inon,1)),squeeze(sum(TwoD.Q(ikm,inon,[4:6]),3)) );
       %using inital profs 
            [ttemp,thtemp,ptemp,lnb2d_vap(ikm,1:length(inon))]=lnb_ice2(T(ikm,inon),P(ikm,inon),pref(:,1),tref(:,1),GridDan(1).Z,squeeze(TwoD.Q(ikm,inon,1)),squeeze(sum(TwoD.Q(ikm,inon,[4:6]),3)) );
   
        end
		if length(inon2)>0
             [ttemp,thtemp,ptemp,lnb2d_tot(ikm,1:length(inon2))]=lnb_ice2(T(ikm,inon2),P(ikm,inon2),pref(:,1),tref(:,1),GridDan(1).Z,squeeze(TwoD.Q(ikm,inon2,1)),squeeze(sum(TwoD.Q(ikm,inon2,[4:6]),3)) );
		end
        
    end    

end

lnb2d_vap=lnb2d_vap/1000+add_ground_height;
lnb2d_tot=lnb2d_tot/1000+add_ground_height;


'done 2d_LNB'

