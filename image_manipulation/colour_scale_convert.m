%converts an image from one colour scale to another from a second image
%notes: some pictures are loaded by imread as MxNx3 arrays where the x3 part is the RGB colour component (and there is no colourmap).
%others are loaded as MxN with a colormap - this is a 3xN array with the x3 part again being RGB. The numbers in the MxN image refer to the 
%row above that represented in the colormap (i.e. to an individual colour). So if number is i then the corresponding colour is in index i+1 of
%the colormap. Numbers of image are in unit format so need to use double command first to perform any arthimetic on them. 

clear c cnew cold

pat_scale='C:\Documents and Settings\Login\My Documents\logbook\vapour_paper\pics\radar_slices\Radar reflectivity (dBZ) 2d plot radar_1_2.png';  
%path for image to use the colour scale from

pat_convert='C:\Documents and Settings\Login\My Documents\logbook\vapour_paper\pics\radar_slices\Fig5b_cpz_24feb04_1753_final.png'; %path for image to be converted

patout='C:\Documents and Settings\Login\My Documents\logbook\vapour_paper\pics\radar_slices\Fig5b_cpz_24feb04_1753_final_col_converted.png'

method='constant_spacing';
%method='look_for_new';


%first image
%[im map]=imread(pat_scale,'gif');
[im,map]=imread(pat_scale);
%image(im(i:i+200,ix-200,ix+200));
ix=820; %x pos of colorbar
i=639; %y pos for bottom of cbar(first colour)
di=44;

c=get_cbar_colours(im,map,i,ix,di,method);
cc=double(c(:,3:end-1,:))/255; %convert from unit8 to double and scale so that max=1
vals=[10:5:70];

%second image
[im2,map2]=imread(pat_convert);
ix=750; %x pos of colorbar
i=1100; %y pos for bottom of cbar(first colour)
di=13;

c2=get_cbar_colours(im2,map2,i,ix,di,method);
cc2=c2(2:33);
vals2=[10:1.5:56.5];

map22=map2;

for i=1:size(vals2,2)
    icol=find(vals2(i)>=vals);
    newcol=cc(1,icol(end),:); %the corresponding colour from the other colourbar
    map22(double(cc2(i))+1,:)=newcol; %replace the colour in the colourmap (need to replace one above the index that would assume)
end

back_col=double(im2(20,10));
map22(back_col+1,:)=[0.7 0.7 0.7];
image(im2);colormap(map22);

back_col=double(im2(1,100));
map22(back_col+1,:)=[1 1 1];
figure;image(im2);colormap(map22);

imwrite(im2,map22,patout,'png');
    
   

%c(33)=0; %black gridlines-put as NaN as don't know values
%values(33)=0; %black