%cube_princ_ax.m
%draws a cube with origin at a corner, the principal axes with the 
%axis of symmetry
clear;
axis ([0,1,0,1,0,1])
%grid on;
box on
view([1 1/2 1/2])%viewpoint in x,y,z cartesian coords
%view(azimuth,elevation)%if desired, maximum is 90 degrees each)
%principal axes
v=1/2; %view best if translate axes by v, other than the symmetry axis
a=line([0,1],[0,1],[0,1],'color', 'r', 'linewidth', 2);%symmetry axis 
%first perpendicular axis and its mirror image
b(1)=line([0+v,1+v],[0+v,-1+v],[0+v,0+v],'LineStyle','-.','color'...
    ,'b','linewidth',2);
line([0+v,-1+v],[0+v,1+v],[0+v,0+v],'LineStyle','-','color','b',...
    'linewidth',2);
%2nd perpendicular axis and its mirror image
b(2)=line([0+v,0+v],[0+v,1+v],[0+v,-1+v],'LineStyle','--',...
    'color','g','linewidth',2);
line([0+v,0+v],[0+v,-1+v],[0+v,1+v],'LineStyle','-','color'...
    ,'g','linewidth',2);
%3rd perpendicular + mirror image, equivalent to 2nd by rotation symmetry
b(3)=line([0+v,1+v],[0+v,0+v],[0+v,-1+v],'LineStyle',':','color','k',...
    'linewidth',2);
line([0+v,-1+v],[0+v,0+v],[0+v,1+v],'LineStyle','-','color',...
    'k','linewidth',2);
title('Unit Cube and Principal Axes','FontSize',13)
xlabel('x','FontSize',14,'Position',[7.6 4.55 3.5]);
ylabel('y','FontSize',14,'Position',[7.6 3.6 3.15]);
zlabel('z','FontSize',14),text(.1,.04,0,'O')
h=legend([a,b],'Princ. Sym. Axis (1,1,1)','2nd Princ. Axis (1,-1,0)',...
    '3rd Princ. Axis (0,1,-1)','Equiv. to 3rd (1,0,-1)',0);
set(h,'FontSize',8,'Position',[0.68 0.83 0.26 0.14])
hold on %let's draw a polygon at the axes intersection with cube edges
x = [1;1/2;0;0;1/2;1;]; y=[0;0;1/2;1;1;1/2;]; z=[1/2;1;1;1/2;0;0];
h=fill3(x,y,z,[0.75 0.75 0.75]);          %draws the hexagon
set(h,'EdgeAlpha',[.3],'FaceAlpha',[0.5]) %hexagon edges, transparent
