%add numbers above each point in timeseries

xnums=xdat(1).x;
ynums=ydat(1).y;
Nnums=num_dat(1).dat;

ynums(round(num_dat(1).dat)==0)='';
xnums(round(num_dat(1).dat)==0)='';
Nnums(round(num_dat(1).dat)==0)='';

clear txt_array
for itxt=1:length(Nnums)
    num_str=num2str(round(Nnums(itxt)));
    txt_array(itxt,1:length(num_str))=num_str;
end

ylims_txt = get(gca,'ylim');
if exist('ioverride_DY') & ioverride_DY==1
    clear ioverride_DY
else
    DY=diff(ylims_txt)/80; %a constant proportion of ylim - useful if not changing ylim
    DY=ynums*0.05; %a percentage of the value (useful when changing ylims)
    DY=0.005; %use a constant value - useful if know final ylim
end

text(double(xnums),double(ynums+DY),txt_array,'horizontalalignment','center','fontweight','bold');