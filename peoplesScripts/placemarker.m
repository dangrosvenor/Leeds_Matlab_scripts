function[ ]=placemarker(fig_h,marker_color,m_size,m_num)
% This function puts marker on the several plots of the same axis
% It also mput markers on legends
% use placemarker
% 
%    [ ]=placemarker(fig_h,marker_color,m_size,m_num)
% where fig_h=figure handle
%      marker_color=color of the marker  
%        m_size=size of the marker 
%        m_num= number of markers to be place; It must be less than plotted data length 
%   Author
% Rajkishore Prasad
% NAIST, Japan
if nargin==0
    fig_h=gcf;
end

figch=get(fig_h,'Children');
L1=length(figch);
ax_h=get(figch(L1),'Children');
L2=length(ax_h);
ind2=1;
for ind1=1:1:L2
    typ=get(ax_h(ind1),'Type');
    if typ=='line'
        grph(ind2)=ax_h(ind1);ind2=ind2+1;
    end
end


L3=length(grph);
for L4=1:1:L3
    
    
    % Marker setting 
    
    colors = strcat('bgrcmyk');
    markers = '+o*.xsdv^><ph'
    if nargin<2
        
        KC =menu('Choose a color for the marker ','Blue','Green','Red','Cyan','Magenta','Yellow',...
            'Black')
        m_color=colors(KC);
        
        MS=menu('Choose a style for the marker ','Plus sign(+)','Circle (o)','Astriscs(*)','Point(.)',....
            'Cross(x)','Square(s)','Diamond(d)','Triangle-down(\/)','Triangle-up(^)','Triangle-right(>)',.....
            'Triangle-left(<)','Five pointed Star(p)','Six pointed Star(h)')
        m_style=markers(MS);
        m_size=input('input marker size>>>')
        m_num= input('input number of marker>>>')%10;% default
    else
        for ind=1:length(marker_color)
            col_indx = strfind(colors,marker_color(ind));
            if ~isempty(col_indx)        %It's a color
                m_color = marker_color(ind);
            else                        %Try a marker instead
                mark_indx = strfind(markers,marker_color(ind));
                if ~isempty(mark_indx)
                    m_style = marker_color(ind);
                end;
            end;
        end;
    end;
    
   % Marker placing on the selected plot
        h=grph(L4);
        m_color=get(h,'Color');
        x1=get(h,'XData');
        y1=get(h,'YData');
        
        Lx1=length(x1);
        Ly1=length(y1);
        if Lx1<m_num
            error('Number of marker is more than data length')
        end
        
        xGap=floor(Lx1/m_num);
        dp=1:xGap:Lx1;
        
        xv=x1(dp);yv=(y1(dp))
        H2 = line(xv,yv,'Color',m_color,'Marker',m_style,'LineStyle','none');
        set(H2,'Markersize',m_size);
    
    
    
    
    
    
    
    
    
    
    
end   %for L4





