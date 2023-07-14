%Glass Water Glass Water Glass Water Glass Water
DMT_dat_1=...
    [0.5 0.6 13.5 11.6 26.5 21.5 39.5 31.8
1.0 0.9 14.0 12.0 27.0 21.9 40.0 32.5
1.5 1.3 14.5 12.5 27.5 22.3 40.5 32.9
2.0 1.8 15.0 12.8 28.0 22.7 41.0 33.3
2.5 2.3 15.5 13.2 28.5 23.1 41.5 33.7
3.0 2.7 16.0 13.6 29.0 23.4 42.0 34.1
3.5 3.1 16.5 13.9 29.5 23.8 42.5 34.5
4.0 3.5 17.0 14.3 30.0 24.4 43.0 34.9
4.5 4.0 17.5 14.7 30.5 24.8 43.5 35.3
5.0 4.5 18.0 14.9 31.0 25.1 44.0 35.7
5.5 5.0 18.5 15.3 31.5 25.5 44.5 36.1
6.0 5.2 19.0 15.7 32.0 25.9 45.0 36.5
6.5 5.5 19.5 16.1 32.5 26.3 45.5 36.9
7.0 5.8 20.0 16.1 33.0 26.7 46.0 37.3
7.5 6.0 20.5 16.5 33.5 27.1 46.5 37.7
8.0 6.2 21.0 17.0 34.0 27.5 47.0 38.2
8.5 6.6 21.5 17.4 34.5 27.9 47.5 38.6
9.0 7.3 22.0 17.8 35.0 28.3 48.0 39.0
9.5 8.1 22.5 18.2 35.5 28.6 48.5 39.4
10.0 8.9 23.0 18.6 36.0 29.0 49.0 39.8
10.5 9.1 23.5 19.0 36.5 29.4 49.5 40.2
11.0 9.5 24.0 19.4 37.0 29.8 50.0 40.6]';
DMT_dat_2=[...
11.5 9.9 24.5 20.0 37.5 30.3
12.0 10.3 25.0 20.4 38.0 30.7
12.5 10.8 25.5 20.8 38.5 31.1
13.0 11.2 26.0 21.2 39.0 31.5]';


DMT_dat=[DMT_dat_1(:); DMT_dat_2(:)];
dmt_glass=DMT_dat(1:2:end);
dmt_water=DMT_dat(2:2:end);
[dmt_glass,I]=sort(dmt_glass);
dmt_water=dmt_water(I);
