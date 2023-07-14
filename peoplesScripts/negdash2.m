function negdash(C,h)
% NEGDASH Make negatively-valued contours dashed.
%       NEGDASH(C,H) uses the output C and H of a contouring 
%       routine such as MATLAB's CONTOUR or Rich P.'s 
%       EXTCONTOUR to reset the LineStyle property of negative
%       contours to '--'. 
%
%       See also CONTOUR, CONTOURC, EXTCONTOUR*, CONTOURSURF*
%       (* by Rich Pawlowicz)

%       For applications where the value of the contour is 
%       repeatedly referenced, you might consider modifying
%       this routine to store the contour value as the UserData
%       property of each line. 

%       Dan Goldner 1/95
%       ----------------------------------------------------

        linehandle=h';
%            x=get(linehandle,'Xdata');startx=x(1);
%            y=get(linehandle,'Ydata');starty=y(1);
%            col=find((C(1,:)==startx))-1;
%            
%            for i=1:length(col)
%                contourvalue=C(1,col(i));
% 
%                if(contourvalue<0),
%                    set(linehandle,'LineStyle',':');
%                end;
%            end
           
           jc=1;
           while jc<size(C,2)
            %ch(1,jc)=str2num(num2str(10^(cbfA(i).c(1,jc)),'%2.2e'));
            
   
                contourvalue = C(1,jc);
                if(contourvalue<0),
                   set(linehandle,'LineStyle',':');
                else
                    set(linehandle,'LineStyle','-');
                end;
                            
            jc=jc+C(2,jc)+1;
        end
           
