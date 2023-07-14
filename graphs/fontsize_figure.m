function fontsize_figure(gcf,gca,fsize,hc)
%fontsize_figure(gcf,gca,fsize)
%changes the fontsize of everything in a figure

%change the legend fontsize - although this removes the marker labels from
%the legend lines...
%set(findobj(gcf,'Type','text'),'FontSize',fsize);

%set(gca,'fontname','times');

%change the x and y label
h_xlabel = get(gca,'XLabel'); 
set(h_xlabel,'FontSize',fsize);
h_ylabel = get(gca,'YLabel'); 
set(h_ylabel,'FontSize',fsize);

%change the ticklabel sizes
set(gca,'FontSize',fsize);

if nargin==4
    set(hc,'Fontsize',fsize);
end