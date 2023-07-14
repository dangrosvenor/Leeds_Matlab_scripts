function [del H0_Houghton]=Houghton_find_thi_streamline_from_mountain_func_version(h,U,gd)


[H_mountain,imax]=max(h);
inds=find(h(1:imax)>1);

[H0_Houghton]=mountain_Houghton_solve_H0_for_h(H_mountain,U,gd);


%del=zeros(length(h));
thi=H0_Houghton*ones([1 length(h)]);
%y_stream=thi;


hx=h(inds);
for i=1:length(inds)
    thi(inds(i))=mountain_Houghton_solve_thi_for_given_h(hx(i),U,H0_Houghton,gd,'windward');
end

hx22=h(imax+1:length(h));
inds=find(hx22>1);
inds0=find(hx22<=1);
hx2=hx22(inds);

L_thi=length(thi);
for i=1:length(inds)
    thi(inds(i)+imax)=mountain_Houghton_solve_thi_for_given_h(hx2(i),U,H0_Houghton,gd,'lee');
end

thi(inds0+imax)=thi(inds(end)+imax);
y_stream = thi + h; %the height of the fluid above z=0
del = y_stream - H0_Houghton;





