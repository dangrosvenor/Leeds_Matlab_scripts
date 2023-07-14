function b=ceil2(a, n)
%function b=ceil2(a, n)
%rounds a to n sigfigs

i0=find(a==0);
in0=find(a~=0);

	aa=abs(a(in0));
	N=log10(aa);
	nn=ceil(abs(N));
	s=sign(a(in0));
	
	btemp=aa./10.^(nn.*sign(N));

    
    altone=find(abs(a(in0))<1);%  & a(in0)>0); %inlcude a>0 as for minus numbers was giving one less sigfig than requested
    nr(altone)=n-1;
    
    altone2=find( abs(a(in0))>1 ); %for a<0 numbers use original n 
    nr(altone2)=n;
    

    
	
	btemp=round2(btemp,nr) .* 10.^(nn.*sign(N));
	
	btemp=btemp.*s;
   
    b(in0)=btemp;
    b(i0)=0;
    
