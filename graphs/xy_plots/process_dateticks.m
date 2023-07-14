function process_dateticks(gca)

tick_text = get(gca,'xticklabel');

for i=1:size(tick_text,1)
    if strfind(tick_text(i,:),'00:00:00')
        str = [tick_text(i,1:6)];
    else
        str = [tick_text(i,13:14)];
    end
    
    for j=1:length(tick_text(i,:))
        tick_text(i,j)=' ';
    end
    pos=length(tick_text(i,:))/2-length(str)/2+1;
    tick_text(i,pos:pos+length(str)-1)=str;
    
end

set(gca,'xticklabel',tick_text);

