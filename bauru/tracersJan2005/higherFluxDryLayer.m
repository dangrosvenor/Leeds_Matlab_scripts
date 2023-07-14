iend=43;
for i=1:length(m)
    m(i).m(:,iend)=sum(m(i).m(:,1:4),2);
end
