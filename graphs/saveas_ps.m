function saveas_ps(fig_handle,filename)
%function saveas_ps(fig_handle,filename) adds the .eps extension
%so don't include this
%Saves the supplied figure handle as a .ps and does the conversion
%of the finsihed file to make the lines in line plots look
%better (using fix_lines script)

filename_ps=[filename '.eps'];
set(fig_handle,'PaperPositionMode','auto');
print(fig_handle,filename_ps,'-depsc');
%fixPSlinestyle(filename);
fix_lines(filename_ps,filename_ps);
