% colour map import for echo top pictures


scrsz=get(0,'ScreenSize');
posit=[9 50 scrsz(3)/5 scrsz(4)/1.2];
figure('Position',posit);

vals=[13 11.5 10 8.5 7 5.5 4 2.5 1];

dire='c:\documents and settings\user\my documents\HIBISCUS\baurufield\radar\echotops\25-28.01\04012500.01C.gif';
[im map]=imread(dire,'gif');
image([11:-1:3]');
colormap(map);

set(gca,'XTicklabel','');
set(gca,'ytick',[1:10]);
set(gca,'yticklabel',vals);
