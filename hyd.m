M = 28*1.67E-27;
k = 1.38E-23;
G = 9.81;

t=35+273.15;
p(1)=1000;
dz=1;

zend=5000;
iend=ceil(zend/dz);

for i=2:iend
    rho=p(i-1)*M/k/t(i-1);
    p(i)=p(i-1)-rho*G*dz;
end


