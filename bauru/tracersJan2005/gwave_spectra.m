clear cq sq qr_pos qr_neg qi_pos qi_neg Epos Eneg E2

iz=132;
tend=1000;
for idir=1:1
    
    for i=1:tend
 
     %   dy=diff(GridDan(idir).Y1(1:2));
	%	wind=TwoD_alltim(idir).W(iz,:,i);
        wind=w2(i,:);
        
        f1=fft(wind);
        cq(i,:)=2*real(f1);
        sq(i,:)=2*imag(f1);
        
    end
    
    cq=cq(:,2:size(cq,2)/2 + 1);
    sq=sq(:,2:size(sq,2)/2 + 1);
    
    for i=1:size(cq,2)
        f1=fft(cq(:,i));
        f2=fft(sq(:,i));
        
        f1=f1(2:size(cq,1)/2);
        f2=f2(2:size(cq,1)/2);
    
        qr_pos(:,i)=0.5 * ( real(f1) - imag(f2) );
        qr_neg(:,i)=0.5 * ( real(f1) + imag(f2) );
        
        qi_pos(:,i)=0.5 * ( real(f2) + imag(f1) );
        qi_neg(:,i)=0.5 * ( real(f2) - imag(f1) );
        
        Epos(:,i) = 2 * (qr_pos(:,i).^2 + qi_pos(:,i).^2);
        Eneg(:,i) = 2 * (qr_neg(:,i).^2 + qi_neg(:,i).^2);
        
    end
    
end

 E=[fliplr(Eneg) Epos];
 
'finished specta'

k=2*pi./[size(Epos,2):-1:1]/1000 ;
w=2*pi./[size(Epos,1):-1:1]/300 ;

w=[1:size(Epos,1)]*2*pi/size(Epos,1)/300;
k=[1:size(Epos,2)]*2*pi/size(Epos,2)/1000;




k=[fliplr(-k) k];

ji=0;
for j=1:2:size(E,1)
    ji=ji+1;
    ii=0;
	for i=1:2:size(E,2)
        ii=ii+1;
        E2(ji,ii)=mean(E(ji,i:i+1),2);
    end
    ii=0;
    for i=1:2:size(E,2)
        ii=ii+1;
        E2(ji+1,ii)=mean(E(ji,i:i+1),2);
    end
    E2(ji,:)=mean(E2(ji:ji+1,:),1);
end


        