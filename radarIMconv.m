dire='c:\documents and settings\user\my documents\HIBISCUS\baurufield\radar\30.01\PPI\2\';
lis=dir(dire);
[a b]=size(lis);
for i=3:a
    fi=strcat(dire,lis(i).name);
    [x y]=size(lis(i).name);
    if y==12
        [im map]=imread(fi,'gif');
        imwrite(im,map,strcat(lis(i).name,'.jpg'),'jpg','quality',100);
    end
end
