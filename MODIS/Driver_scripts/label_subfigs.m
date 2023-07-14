
labs={'a','b','c','d','e','f','g','h','i','j','k','l'};

isub=0;
for i=1:xsub
    for j=1:ysub
        isub=isub+1;
        
        if length(hs)>= isub %sometimes a subplot may be missed out - e.g. if have 5 plots in a 3x2 configuration
            ylims=get(hs{isub},'ylim'); dy = ylims(2)-ylims(1);
            xlims=get(hs{isub},'xlim'); dx = xlims(2)-xlims(1);
            lab_pos_x = xlims(1) - dx/20;
            lab_pos_y = ylims(2) + dy/20;
            axes(hs{isub});
            text(lab_pos_x,lab_pos_y,['(' labs{isub} ')']);
        end
    end
end
