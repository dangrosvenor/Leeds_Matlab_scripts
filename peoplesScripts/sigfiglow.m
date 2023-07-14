function b=sigfigup(a, n)

aa=abs(a);
N=log10(aa);
nn=ceil(abs(N));
s=sign(a);

b=aa/10^(nn*sign(N));

if a<1
    nr=n-1;
else
    nr=n;
end

if s==-1
    b=ceil2(b,nr) * 10^(nn*sign(N));
else
    b=fix2(b,nr) * 10^(nn*sign(N));
end

b=b*s;
