function adlwcgm2,t,p    
    
  
qs=qsatw(t,100.0*p) 

       
rho=100.0*p/(287.04*t)    
g=9.81    
R=287.04    
L=2.5e6    
cp=1004.67    

eps=0.622

;gam in g/kg per metre    
    
;gam in g/m3/m    
    

dqldz=(g*qs/(r*t))*(l*eps/(cp*t)-1.0)/(1+eps*L*L*qs/(r*t*t*cp))
dqldz=rho*dqldz*1000.

;dqldz in g/m3 per m    

return,dqldz    
    

end    
    
