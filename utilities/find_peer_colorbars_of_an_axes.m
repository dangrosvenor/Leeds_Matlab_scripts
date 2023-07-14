function Hc1s = find_peer_colorbars_of_an_axes(Ha1)
    Hf = get(Ha1,'parent');
    Haxs = findobj(Hf,'type','axes');
    IsC=false(1,length(Haxs));
    Hc1s=[];

    for i=1:length(Haxs)
        if isa(handle(Haxs(i)),'scribe.colorbar');
            H=handle(Haxs(i));
            if isequal(double(H.axes),Ha1)
                Hc1s=[Hc1s,Haxs(i)];
            end
        end
    end