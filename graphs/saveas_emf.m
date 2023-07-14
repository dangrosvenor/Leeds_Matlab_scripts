function saveas_ps(fig_handle,filename)
%function saveas_ps(fig_handle,filename)
%saves the supplied figure handle as a .ps and does the conversion
%of the finsihed file to make the lines in line plots look
%better (using fix_lines script)

set(fig_handle,'PaperPositionMode','auto');
%filename_ps = [filename '.eps'];
%print(fig_handle,filename_ps,'-depsc');
%fixPSlinestyle(filename_ps);
%fix_lines(filename_ps,filename_ps);

print(fig_handle,[filename '.emf'],'-dmeta');
%saveas(fig_handle,[filename '.fig'],'fig');


