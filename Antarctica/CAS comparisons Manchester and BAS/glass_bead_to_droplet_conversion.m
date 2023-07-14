
%the following is taken from http://www.ral.ucar.edu/staff/tardif/fog/docs/FM_calib_notes.pdf
%and is for a fog instrument - gives a fairly straight line, but not quite
%Looks to be the same as the DMT one for the CAS as printed in the
%PDF of the calibrations done by DMT.
DAT=[0.5
0.6
14.5
12.5
28.5
23.1
42.5
34.5
1.0
0.9
15.0
12.8
29.0
23.4
43.0
34.9
1.5
1.3
15.5
13.2
29.5
23.8
43.5
35.3
2.0
1.8
16.0
13.6
30.0
24.4
44.0
35.7
2.5
2.3
16.5
13.9
30.5
24.8
44.5
36.1
3.0
2.7
17.0
14.3
31.0
25.1
45.0
36.5
3.5
3.1
17.5
14.7
31.5
25.5
45.5
36.9
4.0
3.5
18.0
14.9
32.0
25.9
46.0
37.3
4.5
4.0
18.5
15.3
32.5
26.3
46.5
37.7
5.0
4.5
19.0
15.7
33.0
26.7
47.0
38.2
5.5
5.0
19.5
16.1
33.5
27.1
47.5
38.6
6.0
5.2
20.0
16.1
34.0
27.5
48.0
39.0
6.5
5.5
20.5
16.5
34.5
27.9
48.5
39.4
7.0
5.8
21.0
17.0
35.0
28.3
49.0
39.8
7.5
6.0
21.5
17.4
35.5
28.6
49.5
40.2
8.0
6.2
22.0
17.8
36.0
29.0
50.0
40.6
8.5
6.6
22.5
18.2
36.5
29.4
9.0
7.3
23.0
18.6
37.0
29.8
9.5
8.1
23.5
19.0
37.5
30.3
10.0
8.9
24.0
19.4
38.0
30.7
10.5
9.1
24.5
20.0
38.5
31.1
11.0
9.5
25.0
20.4
39.0
31.5
11.5
9.9
25.5
20.8
39.5
31.8
12.0
10.3
26.0
21.2
40.0
32.5
12.5
10.8
26.5
21.5
40.5
32.9
13.0
11.2
27.0
21.9
41.0
33.3
13.5
11.6
27.5
22.3
41.5
33.7
14.0
12.0
28.0
22.7
42.0
34.1];

glass_temp=DAT(1:2:end);
water_temp=DAT(2:2:end);

[glass I]=sort(glass_temp);
water = water_temp(I);

iplot_calb=0;
if iplot_calb==1
    figure
    plot(glass,water);
    hold on
    plot([glass(1) glass(end)],[water(1) water(end)],'r-');
    plot(glass,glass,'k');

    xlabel('Glass Diameter (um)');
    ylabel('Water Diameter (um)');

end

%interpolation example
%w=interp1(glass,water,bead)


tols_g=[2.1 8.1 15.5 19.9 20.6 40]; %bead sizes quoted
tols=[0.5 0.8 1.1 1.4 1.4 2.8]; %tolerances (+/-) of beads

bead=60; %required bead size
tol=interp1(tols_g,tols,bead,'linear','extrap');
min_w=interp1(glass,water,bead-tol,'linear','extrap');
max_w=interp1(glass,water,bead+tol,'linear','extrap');
mean_w=interp1(glass,water,bead,'linear','extrap');

fprintf(1,'\nBead = %.1f, tol=%.1f, mean_w=%.1f, min_w=%.1f, max_w=%.1f\n',...
    bead,tol,mean_w,min_w,max_w);





