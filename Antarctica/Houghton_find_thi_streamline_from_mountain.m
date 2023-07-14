clear thi
H_mountain=1500;

gd=0.7;
U=8;

[H0_Houghton]=mountain_Houghton_solve_H0_for_h(H_mountain,U,gd);

hx=[1:1:H_mountain];

for i=1:length(hx)
    thi(i)=mountain_Houghton_solve_thi_for_given_h(hx(i),U,H0_Houghton,gd,'windward');
end

hx2=[H_mountain-1:-1:1];

L_thi=length(thi);
for i=1:length(hx2)
    thi(i+L_thi)=mountain_Houghton_solve_thi_for_given_h(hx2(i),U,H0_Houghton,gd,'lee');
end

y_stream = thi + [hx hx2]; %the height of the fluid above z=0



