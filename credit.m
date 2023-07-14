x2(1)=1000;
p=1.01093;
m=0.98;
set=30;
n=48;

for i=2:n+1
    
    x2(i)=(x2(i-1)-set)*p;
    int2(i)=(x2(i-1)-set)*(p-1);
end

b=(1-(1/p));
a=(p*m);
sum=0;
sum2=0;

fprintf('%d %f %f %f %f\n',1,x2(1),0,x2(1),int2(1));
for i=2:n+1
    
    x1=a^(i-1)*x2(1);
    int1=x1*b;
    sum=sum+int1;
    sum2=sum2+int2(i);
    
    
    fprintf('%d %f %f %f %f %f %f\n',i,x1,int1,sum,x2(i),int2(i),sum2);
    
end

%fprintf('%d %f %f %f %f\n',0,0,sum,0,sum2);
    