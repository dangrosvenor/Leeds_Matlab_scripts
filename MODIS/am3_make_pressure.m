function am3_phalf=am3_make_pressure(am3_ps,am3_pk,am3_bk)
%normal order of other 3D variables is (time,level,lat,lon)

%make 3D arrays
s_horiz = size(am3_ps);
s_vert = length(am3_pk);

if length(s_horiz)==3
    ps = repmat(am3_ps,[1 1 1 s_vert]);
    pk = repmat(am3_pk,[1 s_horiz(1) s_horiz(2) s_horiz(3)]);
    bk = repmat(am3_bk,[1 s_horiz(1) s_horiz(2) s_horiz(3)]);

    ps=permute(ps,[1 4 2 3]);
    pk=permute(pk,[2 1 3 4]);
    bk=permute(bk,[2 1 3 4]);

else
    ps = repmat(am3_ps,[1 1 s_vert]);
    pk = repmat(am3_pk,[1 s_horiz(1) s_horiz(2)]);
    bk = repmat(am3_bk,[1 s_horiz(1) s_horiz(2)]);

    ps=permute(ps,[3 1 2]);
    %pk=permute(pk,[2 1 3 4]);
    %bk=permute(bk,[2 1 3 4]);
end

am3_phalf = pk + ps.*bk;