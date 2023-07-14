%x_grid = ([1:size(lat2d(1).var,2)]-1) * distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(1,2),lon2d.var(1,2));
%y_grid = ([1:size(lat2d(1).var,1)]-1) * distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(2,1),lon2d.var(2,1));

x_grid=x_grid_orig;
y_grid=zz(1).z;


H=max(y_grid); %draw streamlines up to this height (km)draw_streamlines
H=3.2; %draw streamlines up to this height (km)draw_streamlines

ix=findheight(x_grid,1000);
iy=findheight_nearest(y_grid,H);    


%udat=u_310_2;
%vdat=v_310_2;

%udat=wrf_zint_u;
%vdat=wrf_zint_v;

udat=u_quiver;
vdat=v_quiver;

[X,Y] = meshgrid(x_grid,y_grid);


dx=1;
dy=1;
recalc=1;


switch recalc
    case 1
        clear STREAM
        i=2:dy:iy;
        j=iy:-dy:1;
        count=0;
        %       for i=1:dx:ix
        for js=1:1  %length(j)
           count=count+1;
           STREAM{count}=stream2(X,Y,udat,vdat,x_grid(1),y_grid(j(js)));
        end
        %        end
        
%ih=length(i);
ih=findheight(y_grid,2);
        for js=1:ih
            count=count+1;
%            STREAM{count}=stream2(X,Y,-udat,-vdat,x_grid(end),y_grid(js));
        end

end
%
% x=[];
% y=[];

plot_type='lines';
%plot_type='patch';
switch plot_type
    case 'lines'
%         switch ih_wrf
%             case 10
%                 col='r';
%             case 32
%                 col='b';
%         end

        col='b';
        
        for i=1:length(STREAM)
            plot(STREAM{i}{1}(:,1),STREAM{i}{1}(:,2),col,'linewidth',1.3);
            % %    x=[x; STREAM{i}{1}(:,1)];
            % %    y=[y; STREAM{i}{1}(:,2)];
        end

    case 'patch'
        switch ih_wrf
            case 10
                col=[0.5 0 0];
            case 32
                col=[0 0 0.5];
        end
        patx=[STREAM{1}{1}(:,1); flipud(STREAM{end}{1}(:,1)); 0];
        paty=[STREAM{1}{1}(:,2); flipud(STREAM{end}{1}(:,2)); 0];
        h=patch(patx,paty,col);
        alpha(h,0.25);
end


