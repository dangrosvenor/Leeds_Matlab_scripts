%Make the maps of LWP at night, coarse grain to AMSRE resolution and then
%calculate the cloud fraction based on LWP>20

%This makes the maps and coarsens to 0.25 degrees
UM_maps_LWP_NIGHT_time_loop_RUN_v1

for i=1:length(UM_time_out)
    
    icf = find(UM_time_out{i}.datUM{1}>=20);
    CF(i) = length(icf)./length(UM_time_out{i}.datUM{1}(:));
    
end