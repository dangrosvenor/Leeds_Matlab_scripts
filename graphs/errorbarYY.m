function [hline,bar_end,bar_beg]=errorbarYY(horiz_vert,x,y,errorb,parent,color,marker,LW,end_cap_fraction)
%function[hline]=errorbarYY(x,y,errorb,parent,color,marker,LW,end_cap_size)
% LW is the linewidth
% end_cap_fraction is the fraction of the whole span of the plot ('xlim') that the
% end cap should be 
%set parent to the axis to be plotted on
% example - herr=errorbarYY('vert',1847.5,me_t_PI,std_t_PI*2,gca,'r','o',2,0.01);
%Written by Dan McCoy and Dan  Grosvenor

switch horiz_vert
    case 'horiz'
        lim_size = get(parent,'ylim');
%        lim_size = get(parent,'xlim');
    case {'vert filled bar','vert'}        
        %lim_size = [0 spanX];
        %lim_size = [0 spanX];
        lim_size = get(parent,'xlim');
end

end_cap_size = diff(lim_size)*end_cap_fraction;


%LW=2;
%end_cap_size=5;
%errorb
%hline=line(x,y,'Color',color,'Parent',parent,'LineWidth',LW,'Marker',marker);
%% WHY THE **(#@)( CAN'T MATLAB DO THIS

switch horiz_vert
    case 'horiz'
        for i=1:length(x)

            bar_end=x(i)-errorb(i);
            bar_beg=x(i)+errorb(i);

            hline((i-1)*3+1)=line([bar_end
                bar_beg],[y(i) y(i)],'Color',color,'Parent',parent,'LineWidth',LW,'LineStyle','-');
            hline((i-1)*3+2)=line([bar_end
                bar_end],[y(i)-end_cap_size y(i)+end_cap_size],'Color',color,'Parent',parent,'LineWidth',LW,'LineStyle','-');
            hline((i-1)*3+3)=line([bar_beg
                bar_beg],[y(i)-end_cap_size y(i)+end_cap_size],'Color',color,'Parent',parent,'LineWidth',LW,'LineStyle','-');
        end


    case 'vert'
        for i=1:length(x)

            bar_end=y(i)-errorb(i);
            bar_beg=y(i)+errorb(i);

            hline((i-1)*3+1)=line([x(i) x(i)],[bar_end
                bar_beg],'Color',color,'Parent',parent,'LineWidth',LW,'LineStyle','-');
            hline((i-1)*3+2)=line([x(i)-end_cap_size x(i)+end_cap_size],[bar_end
                bar_end],'Color',color,'Parent',parent,'LineWidth',LW,'LineStyle','-');
            hline((i-1)*3+3)=line([x(i)-end_cap_size x(i)+end_cap_size],[bar_beg
                bar_beg],'Color',color,'Parent',parent,'LineWidth',LW,'LineStyle','-');
        end
        
    case 'vert filled bar'
        for i=1:length(x)

            bar_end=y(i)-errorb(i);
            bar_beg=y(i)+errorb(i);              

           X = [x(i)-end_cap_size x(i)-end_cap_size x(i)+end_cap_size x(i)+end_cap_size];
           Y = [bar_beg bar_end bar_end bar_beg];
           
           hline = fill(X,Y,color,'linestyle','none');

        end

end

'';





