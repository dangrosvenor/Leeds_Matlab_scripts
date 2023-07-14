%echotops
val_tops=[5:70];

tops=zeros(size(zh(:,:,1)));

for itop=1:length(val_tops)
    
	for iz=1:length(zar)
        i10=find( zh(:,:,iz)>=val_tops(itop) );
        tops(i10)=zar(iz);
        ntops(iz,itop,iarm_radar)=length(i10); %number of points exceeding x dbz
	end

end

