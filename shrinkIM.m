%imports several pictures and makes them smaller

patt='c:\documents and settings\user\my documents\HIBISCUS\baurufield\photos\recup_nilocube\newspaper\';
pout='c:\documents and settings\user\my documents\HIBISCUS\baurufield\photos\matjpgs\';

lis=dir(patt);
[aeta beta]=size(lis);

for ieta=3:aeta
    pat=strcat(patt,lis(ieta).name);
    [siza sizb]=size(lis(ieta).name);
    if sizb==12
        if strcmp('JPG',lis(ieta).name(10:12))==1
            im=imread(pat,'jpg');
            pout2=strcat(pout,lis(ieta).name);
            imwrite(im,pout2,'jpg','quality',30);
        end
    end
end
