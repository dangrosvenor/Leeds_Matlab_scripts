function me=nzeromean(x,dim)

x=permute(x,[dim 3-dim]); %reoder matrix so dim dimension is the first

b=size(x);
ndims=length(b);

for i=1:b(1)
    in0=find(x(i,:)~=0);
    me(i)=mean(x(i,in0),2);
end

