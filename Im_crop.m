%Im_crop

pat='c:/documents and settings/g/my documents/my pictures/me_and_lou_glastonbury/';
pat2='c:/documents and settings/g/my documents/my pictures/glasto_cropped/';

lis=dir(pat);


%Lou
Lx=600*0.85;
Ly=590*0.85;

ix=815;
iy=525;

coords(3).c=[num2str(iy) ':' num2str(iy+Ly) ',' num2str(ix) ':' num2str(ix+Lx) ',:'];

%Dan
Lx=600*1.15;
Ly=590*1.15;

ix=775;
iy=355;

coords(4).c=[num2str(iy) ':' num2str(iy+Ly) ',' num2str(ix) ':' num2str(ix+Lx) ',:'];

%Dave
Lx=700;
Ly=500;

Lx=600*1.05;
Ly=590*1.05;

ix=765;
iy=425;

coords(6).c=[num2str(iy) ':' num2str(iy+Ly) ',' num2str(ix) ':' num2str(ix+Lx) ',:'];


ims=[6];   %6];
for i=ims   %3:length(lis)
    
   % im(i).i=imread(strcat(pat,lis(i+2).name),'jpg');
    figure
    comm=['image(im(i).i(' coords(i).c '))'];
    eval(comm);
    outname=[pat2 'out_' num2str(i) '.jpg'];  
    eval(['im_out=im(i).i(' coords(i).c ');']);
end





imwrite(im_out,outname,'jpg');     %,'quality',25);
