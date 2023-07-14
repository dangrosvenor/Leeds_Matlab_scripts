function b6=LiuB6f(dis)

%dis=[0.01:0.01:0.5];
alx=dis.^-2 - 1;

b6=(gamma(6+alx+1)./gamma(1+alx)).^(1/6) .* (gamma(3+alx+1)./gamma(1+alx)).^(-1/3);

%figure;
%plot(dis,b6.^6);