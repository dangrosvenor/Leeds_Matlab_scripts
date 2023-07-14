function saveas_ps_fig_emf(fig_handle,filename,print_filename,tag,nofixlines_ps)
%function saveas_ps(fig_handle,filename,print_filename,tag,nofixlines_ps)
%saves the supplied figure handle as a .ps and does the conversion
%of the finsihed file to make the lines in line plots look
%better (using fix_lines script)
%If don't want to fix the lines (makes some color plots go wrong) then
%input nofixlines_ps as 1 (otherwise zero, or leave out).

if nargin<3
    print_filename=1;
end

if nargin>=4
    filename=[filename tag];
end

if nargin<5
    nofixlines_ps=0;
end

set(fig_handle,'PaperPositionMode','auto');

ylims=get(gca,'ylim');
xlims=get(gca,'xlim');
ydir=get(gca,'ydir');
if strcmp(ydir,'reverse')
    ylims_save=ylims(1);
    ylims(1)=ylims(2);
    ylims(2)=ylims_save;
end
ypos = ( ylims(2) - ylims(1) )*0.07 + ylims(2);
xpos = xlims(1) - ( xlims(2) - xlims(1) )*0.1;



if print_filename==1
    text(xpos,ypos,filename,'fontsize',8);
end




filename_ps = [filename '.eps'];
%print(fig_handle,filename_ps,'-depsc','-tiff','-r150');
print(fig_handle,filename_ps,'-depsc2');

if nofixlines_ps==0
    fixPSlinestyle(filename_ps);
end
%fix_lines(filename_ps,filename_ps);

print(fig_handle,[filename '.emf'],'-dmeta');
print(fig_handle,[filename '.png'],'-r600','-dpng');
saveas(fig_handle,[filename '.fig'],'fig');

if print_filename==1
   unplot %remove the text
end



