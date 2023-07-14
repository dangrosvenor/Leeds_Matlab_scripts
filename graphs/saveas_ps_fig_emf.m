function [filename]=saveas_ps_fig_emf(fig_handle,filename,tag,itext,idate,ititle_nice,title_nice,save_notes_filepath,iplot_eps,optional)
% function [filename]=saveas_ps_fig_emf(fig_handle,filename,tag,itext,idate,ititle_nice,title_nice,save_notes_filepath,iplot_eps,optional)
% Or just run using saveas_ps_fig_emf(fig_handle,filename)
% saves the supplied figure handle as a .ps and does the conversion
% of the finsihed file to make the lines in line plots look
% better (using fix_lines script)
% itext tells it whether to put the full path filename at the top of image
% (for location purposes)
% idate is whether to put a unique date/time on the end of filename to avoid
% overwriting
% optional is a structure

if exist('optional')==1
    %Convert all of the variable names in the input structure to actual names
    %for ease of use
    name_struc='optional'; %The name of the structure
    names = eval(['fieldnames(' name_struc ');']);
    for i=1:length(names)
        eval_str = [names{i} ' = ' name_struc '.' names{i} ';'];
        eval(eval_str);
    end
end

%this is better for eps as is produces scalable fonts, etc. But there is
%the problem of criss-crossing lines - although this is possibly the fault
%of anti-aliasing in the viewer (Adobe, etc). Not sure if it can be fixed
%or not so that viewers by default don't show the anti-aliasing??
set(fig_handle,'renderer','painters');
%also may take a while with complicated images. Perhaps openGL (the default
%might be better in that case).

if ~exist('iplot_eps')
    iplot_eps=1; %flag to plot eps (and pdf) output
end
if ~exist('isavefig')
    isavefig=1; %Flag for whether to save as .fig
end
if ~exist('iplot_png')
    iplot_png=0; %Flag for whether to save as .png using Matlab method - but this does not work as well as the one via the eps (small font - seems related to resolution requested)
end
if ~exist('iplot_jpg')
    iplot_jpg=1; %Flag for whether to save as .fig
end

if ~exist('tag')
    tag='';
end

if ~exist('itext')
    itext=1;
end

if ~exist('idate')
    idate=1;
end

if ~exist('isave_data')
    isave_data=0;
end


if exist('ititle_nice')
    if ititle_nice==1
        h=get(gca,'title');
        title_save = get(h,'string');
        title(title_nice);
    end
else
    ititle_nice = 0;
end

if exist('save_notes_filepath')
    isave_notes=1;      
else
    isave_notes=0;
end

set(fig_handle,'PaperPositionMode','auto');

ylims=get(gca,'ylim');
xlims=get(gca,'xlim');
ypos = ( ylims(2) - ylims(1) )*0.15 + ylims(2);
xpos = xlims(1) - ( xlims(2) - xlims(1) )*0.1;

if itext==1
    text(xpos,ypos,filename,'fontsize',8);
end

LT=length(tag);
Nmax=200;
if length(filename)+LT>Nmax
    filename(Nmax-LT:end)='';
end
if idate==1
    datestr_now = datestr(now,30);
    filename=[filename '_' datestr_now];
end
    
filename = [filename tag];

%print(fig_handle,filename_ps,'-depsc','-tiff','-r150');
filename=remove_character(filename,'*','');
filename=remove_character(filename,'\','');
filename=remove_character(filename,' ','_'); %replace all spaces with underscores - latex can handle single
filename=remove_character(filename,'<','.LT.'); %replace all spaces with underscores - latex can handle single
filename=remove_character(filename,'>','.GT.'); %replace all spaces with underscores - latex can handle single
filename=remove_character(filename,'%','pct'); %replace all spaces with underscores - latex can handle single spaces, but not double - may as well remove all spaces for now
filename=remove_character(filename,':',''); %
filename=remove_character(filename,'{',''); %
filename=remove_character(filename,'}',''); %


%keep this line - used by other parts
filename_ps = [filename '.eps'];
istr = findstr(filename_ps,'/');

if iplot_eps==1
    print(fig_handle,filename_ps,'-depsc2');

    %fixPSlinestyle(filename_ps);   %fix_lines(filename_ps,filename_ps);
    eval(['!epstopdf ' filename_ps(1:istr(end)) '''' filename_ps(istr(end)+1:end) '''']);

    %the png output produces a very small font, so this seems to work better:-
    eval(['!convert +antialias -density 350 ' filename_ps(1:istr(end)) '''' filename_ps(istr(end)+1:end) ''' ' filename_ps(1:istr(end)) '''' filename_ps(istr(end)+1:end-4) '.png'';']);
%    eval(['!convert +antialias -density 600 ' filename_ps(1:istr(end)) '''' filename_ps(istr(end)+1:end) ''' ' filename_ps(1:istr(end)) '''' filename_ps(istr(end)+1:end-4) '.png'';']);    
    % need -trim option? Prob not.
end

if iplot_png==1
    print(fig_handle,[filename '.png'],'-r350','-dpng');
end

%print(fig_handle,[filename '.svg'],'-dsvg');%Doesn't work on Olympus
    %at least.

if iplot_jpg==1
    %print(fig_handle,[filename '.emf'],'-dmeta'); %Doesn't work on Olympus
    %at least.
    print(fig_handle,[filename '.jpg'],'-r300','-djpeg');
%    print(fig_handle,[filename '.jpg'],'-r800','-djpeg');
end

if isavefig==1
    saveas(fig_handle,[filename '.fig'],'fig');
end




if itext==1
    unplot %remove the text
end

if ititle_nice==1
    title(title_save);
end

if isave_notes==1
    eval(['!cp ' save_notes_filepath ' ' filename '.txt']);
%    fid=fopen([filename '.txt'],'wt');
%    fprintf(fid,'%s',save_notes);
%    fclose(fid);
end

if isave_data==1
    save([filename '.mat'],'plot_data','plot_data_restricted','lon_data','lat_data','lon_data_edges','lat_data_edges','-V7.3');
end

