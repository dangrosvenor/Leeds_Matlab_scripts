function [PRESS,RHO] = Press(TH,THREF,PSF,YI,ZI)


% THETA ON MESH
size(TH);
jmax=ans(1);
size(TH);
kmax=ans(2);
for j=1:jmax
    %for k=1:kmax
        THETA(j,1:kmax)=TH(j,1:kmax)+THREF(1:kmax);
        %end
end
PRESS=zeros(jmax,kmax);
    % DEFINE SOME CONSTANTS
    K = 0.286;
    MPrime = 28*1.67E-27;
    k_boltz = 1.4E-23;
    GRAV = 9.8;
    
PRESS(:,1) = PSF;
for k =1:(kmax-1)
    for j=1:(jmax)      
        m = (THETA(j,k+1)'-THETA(j,k)')./(ZI(k+1,j)-ZI(k,j));
        c = THETA(j,k)'-((THETA(j,k+1)'-THETA(j,k))./(ZI(k+1,j)-ZI(k,j))).*ZI(k,j);
        if(m == 0)
            PRESS(j,k+1) =PRESS(j,k);
        else
            var=log((m.*ZI(k+1,j)'+c)./(m.*ZI(k,j)'+c));
            var = var.*MPrime.*GRAV.*(100000.^K)./(k_boltz.*m);
            PRESS(j,k+1) = ((PRESS(j,k).^K)./K) - var;
            PRESS(j,k+1) = (PRESS(j,k+1).*K).^(1./K);
            PRESS(j,k+1) = real(PRESS(j,k+1));      
        end
    end
end
clear var;
% NOW TO WORK OUT THE DENSITY OF THE AIR
for k = 1:kmax
    %for j = 1:(jmax)
        RHO(1:jmax,k) = real(PRESS(1:jmax,k)*MPrime./(real(THETA(1:jmax,k)).*k_boltz.*((PRESS(1:jmax,k)./100000.).^K)));
        %end
end