patt='c:/cygwin/home/user/lesexe/les2.3/';


lis=dir(patt);
[aeta beta]=size(lis);

for ieta=3:aeta
    pat=strcat(patt,lis(ieta).name);
    [siza sizb]=size(lis(ieta).name);
    
    if findstr(pat,'.f')
        fid=fopen(pat,'rb');
        status=0;
        count=0;
    while status==0 & count<350
        count=count+1;
        
        status=fseek(fid,1,'cof');
        fseek(fid,-1,'cof');
        if status==0
            %lis(ieta).name
            li=fgetl(fid);
            if(findstr(li,'IUSETHP'))
                lis(ieta).name
                break
            end
    
        end
        
    end
    
    end
end
    
    