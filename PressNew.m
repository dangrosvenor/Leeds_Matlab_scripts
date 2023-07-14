function [PRESS,RHO,TEMP,QSAT] = PressNew(GridPress,TwoDPress)

% THETA ON MESH
[YI2 ZI2]=meshgrid(GridPress.Y1,GridPress.Z);
THETA=TwoDPress.TH2;
size(THETA);
dom=ans(2);
dom2=ans(1);
%dom
%dom2

% DEFINE SOME CONSTANTS
    K = 0.286;
    MPrime = 28*1.67E-27;
    k = 1.38E-23;
    GRAV = 9.8;
    

if size(TwoDPress.PP,1)==1

PRESS=zeros(dom2,dom);
RHO=zeros(dom2,dom);

    
   
    
PRESS(1,:) = 94100.;
fprintf(1,'surface pressure = %f  - check for particular sounding',PRESS(1,1));
for i=1:(dom2-1)
    for j=1:dom      
        m = (THETA(i+1,j)-THETA(i,j))/(ZI2(i+1,j)-ZI2(i,j));
        c = THETA(i,j)-((THETA(i+1,j)-THETA(i,j))/(ZI2(i+1,j)-ZI2(i,j)))*ZI2(i,j);
        if(m == 0.)
            PRESS(i+1,j) =PRESS(i,j);
        else
            var=log(abs(m.*ZI2(i+1,j)+c)./(m.*ZI2(i,j)+c));
            var = var.*MPrime.*GRAV.*(100000.^K)./(k.*m);
            PRESS(i+1,j) = ((PRESS(i,j).^K)./K) - var;
            PRESS(i+1,j) = (PRESS(i+1,j).*K).^(1./K);
            PRESS(i+1,j) = real(PRESS(i+1,j));      
        end
   end
end
clear var;
% NOW TO WORK OUT THE DENSITY AND TEMP AND QSAT OF THE AIR
for i = 1:dom2
    for j = 1:dom
        RHO(i,j) = (PRESS(i,j)*MPrime./((THETA(i,j))*k*((PRESS(i,j)./100000.).^K)));
        TEMP(i,j) = THETA(i,j)*(PRESS(i,j)/100000)^K;
        QSAT(i,j) = 3.8/( (PRESS(i,j)/100) *exp(-17.27*(TEMP(i,j)-273.15)/(TEMP(i,j)-35.86)) -6.109 );
        if i==1 
            QSAT(i,j)=0;
        end;
    end
end

else
    
    PRESS=TwoDPress.PP;
    TEMP = THETA.*(PRESS/100000).^K;
    RHO = PRESS.*MPrime./(k*TEMP);
    QSAT = 3.8./( (PRESS/100) .*exp(-17.27*(TEMP-273.15)./(TEMP-35.86)) -6.109 );
    
end


    
    
    