
for ibox=1:length(xbox)

x=xbox(ibox).x(1);
y=ybox(ibox).y(1);

w = xbox(ibox).x(2) - x;
h = ybox(ibox).y(2) - y;

X=[x x+w x+w x];
Y=[y y y+h y+h];



    switch box_type
        case 'fill'
%                fill(X,Y,colbox(ibox).c,'linestyle','none','facealpha',0.3);
                H_fill=fill(X,Y,colbox(ibox).c,'linestyle','none');     
        case 'rectangle'
            %        rectangle('Position',[x/24,y,w/24,h],'LineWidth',2,'LineStyle','--');
            line(X([1 4]),Y([1 4]),'LineWidth',2,'LineStyle','--','color',colbox(ibox).c);
            line(X([2 3]),Y([2 3]),'LineWidth',2,'LineStyle','--','color',colbox(ibox).c);
            line(X([1 2]),Y([3 3]),'LineWidth',2,'LineStyle','--','color',colbox(ibox).c);
            line(X([1 2]),Y([2 2]),'LineWidth',2,'LineStyle','--','color',colbox(ibox).c);
    end
    
    
end