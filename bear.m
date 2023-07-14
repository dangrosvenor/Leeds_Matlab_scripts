function an=bear(u,v)
clear pi

an=abs(atan(u./v).*180./pi);

for i=1:size(u,1)
    
if (u(i)>0 & v(i)<0)
    an(i)=an(i)+90;
elseif (u(i)<0 & v(i)<0)
    an(i)=an(i)+180
elseif (u(i)<0 & v(i)>0)
    an(i)=an(i)+270
end

end