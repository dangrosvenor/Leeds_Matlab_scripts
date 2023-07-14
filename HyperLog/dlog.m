function d=dlog(x,xc)

i=find(x>=xc);
d(i)=log10(x(i));

i=find(x<-xc);
d(i)=2*log10(xc)-log10(-x(i))-1;

i=find(x<xc & x>=-xc);
a=log10(xc);
b=a-1;
d(i)=a+(x(i)-xc)/(-2*xc) * (b-a);

d=reshape(d,size(x));

% 
% for i=1:length(x)
%     
% 	if x(i)>=xc
%         d(i)=log10(x(i));
% 	elseif x(i)<-xc
% 		d(i)=2*log10(xc)-log10(-x(i))-1;
% 	else
%         a=log10(xc);
%         b=a-1;
%         d(i)=a+(x(i)-xc)/(-2*xc) * (b-a);
% 	end
% 
% end

