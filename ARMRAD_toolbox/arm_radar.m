function varargout = arm_radar(varargin)
% ARM_RADAR Application M-file for arm_radar.fig
%    FIG = ARM_RADAR launch arm_radar GUI.
%    ARM_RADAR('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 29-Apr-2006 14:41:45

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

    set(handles.figure1,'renderer','opengl');
    ncst=get(handles.figure1,'userdata');%load([matlabroot,'/toolbox/ARMRAD_toolbox/TiwisLL']);
    axes(handles.axes1);
    plot(ncst(find(ncst(:,2)<0 |isnan(ncst(:,2))),1),ncst(find(ncst(:,2)<0 |isnan(ncst(:,2))),2));
    axis tight;grid on;
    %set(handles.axes1,'zlim',[0 20],'view',[45 45]);hold on;
    set(handles.axes1,'zlim',[0 20]);hold on;
    % Set to be at the centre of the screen
    centreGUI(handles.figure1);

    if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		disp(lasterr);
	end

end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.



% --------------------------------------------------------------------
function varargout = figure1_ResizeFcn(h, eventdata, handles, varargin)
% Stub for ResizeFcn of the figure handles.figure1.
disp('figure1 ResizeFcn not implemented yet.')


% --------------------------------------------------------------------
function varargout = edit1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit1.
disp('edit1 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = pushbutton1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton1.
disp('pushbutton1 Callback not implemented yet.')

fp=get(handles.edit1,'string');
try
    [FILENAMES, PATHNAME] = uigetfiles('*.ascii','Pick arm radar files',fp);
    if(~length(FILENAMES))
        error('no files');
    end
catch
    [FILENAMES, PATHNAME] = uigetfiles('*.ascii','Pick arm radar files','\');
end

if(~length(FILENAMES))
    return;
end

UD_push1.filenames = FILENAMES;
UD_push1.pathname  = PATHNAME;

set(h,'userdata',UD_push1);

% --------------------------------------------------------------------
function varargout = pushbutton2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton2.
disp('pushbutton2 Callback not implemented yet.')

UD_push1 = get(handles.pushbutton1,'userdata');

if(~isfield(UD_push1,'filenames'))
    uiwait(errordlg('Pick files first','error','modal'));
    return;
end

% Import the data
h1 = waitbar(0,'Please wait...');
for i=1:length(UD_push1.filenames)
    [xar,yar,zar,zh,hclass,slatr,slongr,adate,ahhmm,latitude,longitude,radClass]=...
        radarClassARMPaul([UD_push1.pathname,UD_push1.filenames{i}]);
    UD_push2(i).xar = xar;
    UD_push2(i).yar = yar;
    UD_push2(i).zar = zar;
    UD_push2(i).zh  = zh;
    UD_push2(i).hclass = hclass;
    UD_push2(i).slatr = slatr;
    UD_push2(i).slongr = slongr;
    UD_push2(i).adate = adate;
    UD_push2(i).ahhmm = ahhmm;
    UD_push2(i).latitude = latitude;
    UD_push2(i).longitude = longitude;
    UD_push2(i).radClass = radClass;
    waitbar(i/length(UD_push1.filenames),h1);
end
close(h1);
assignin('base','RADDAT',UD_push2);
set(h,'userdata',UD_push2);


% --------------------------------------------------------------------
function varargout = pushbutton3_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton3.
disp('pushbutton3 Callback not implemented yet.')

set(handles.pushbutton4,'userdata',0);
%     load([matlabroot,'/toolbox/ARMRAD_toolbox/TiwisLL']);
%     axes(handles.axes1);
%     plot(ncst(find(ncst(:,2)<0 |isnan(ncst(:,2))),1),ncst(find(ncst(:,2)<0 |isnan(ncst(:,2))),2));
%     axis tight;grid on;
%     set(handles.axes1,'zlim',[0 20],'view',[45 45]);hold on;
% 
% Loop
UD_push2 = get(handles.pushbutton2,'userdata');
if(~length(UD_push2))
    uiwait(errordlg('Import files first','error','modal'));
    return;
end
valSel=get(handles.popupmenu1,'value');

while(1)
    stopFlag=get(handles.pushbutton4,'userdata');
    if(stopFlag ==1)
        break;
    end
for i=1:length(UD_push2)
    stopFlag=get(handles.pushbutton4,'userdata');
    if(stopFlag ==1)
        break;
    end
    UD_push2(i).zh(find(UD_push2(i).zh(:)<-90 )) = NaN;

    % Where do you want to start slice?
%     hFig=figure('visible','off');
%     xCenStart = 0.;
%     yCenStart = 0.;
%     zCenStart = 10.;
%     hsf = surf(UD_push2(i).xar,...
%         yCenStart+linspace(0.,sqrt((max(UD_push2(i).zar).^2)),length(UD_push2(i).yar)), ...
%         zCenStart.*ones(length(UD_push2(i).yar),length(UD_push2(i).yar)));
%     
%     rotate(hsf,[1 0 0],0.01,[0 yCenStart 0]);
%     %rotate(hsf,[0 0 1],45,[0 0.5 0]);
%     xd = get(hsf,'XData');
%     yd = get(hsf,'YData');
%     zd = get(hsf,'ZData');

    xd=UD_push2(i).xar;
    yd=UD_push2(i).yar;
    [xd,yd]=meshgrid(xd,yd);
    zd=6.*ones(size(xd));
    % Convert to lat, long, z
	latd = UD_push2(i).slatr + yd./1.852./60;
	Re=6378.1e3;
	longd = UD_push2(i).slongr + 1./(2.*pi.*Re.*cos(latd.*pi./180)./360) .*xd.*1000;
    %delete(hsf);

    
    axes(handles.axes1);
    switch valSel
    case 1
%         hSlice=slice(UD_push2(i).longitude,UD_push2(i).latitude,UD_push2(i).zar,...
%             UD_push2(i).zh(:,:,:),longd,latd,zd);shading flat;caxis([-10 60]);
%         hc=colorbar;set(hc,'visible','on');
%         set(hc,'yaxislocation','right');
        hSlice=slice(UD_push2(i).longitude,UD_push2(i).latitude,UD_push2(i).zar,...
            UD_push2(i).zh(:,:,:),[100:200],[-11.4],[5]);shading flat;caxis([-10 60]);
        hc=colorbar;set(hc,'visible','on');
        set(hc,'yaxislocation','right');
    case 2
%         hSlice=slice(UD_push2(i).longitude,UD_push2(i).latitude,UD_push2(i).zar,...
%             UD_push2(i).hclass(:,:,:),longd,latd,zd);shading flat;caxis([0 10]);
        hSlice=slice(UD_push2(i).longitude,UD_push2(i).latitude,UD_push2(i).zar,...
            UD_push2(i).hclass(:,:,:),131,-11.66,12);shading flat;caxis([0 10]);
        hc=colorbar;set(hc,'visible','on');
        set(hc,'ytick',[0:1:10]);
        set(hc,'yticklabel',UD_push2(i).radClass,'yaxislocation','right','fontsize',8);
    end
%     hSlice=slice(UD_push2(i).longitude,UD_push2(i).latitude,UD_push2(i).zar,...
%         UD_push2(i).zh(:,:,:),xd,yd,zd);shading flat;caxis([-10 60]);
%        slice(UD_push2(i).xar,UD_push2(i).yar,UD_push2(i).zar,...
%            UD_push2(i).hclass(:,:,:),xd,yd,zd);shading flat;caxis([0 10]);

    title([UD_push2(i).adate,' ',UD_push2(i).ahhmm]);
%    print(gcf, '-djpeg',['picture',num2str(i),'.jpg']);
    pause(0.5);
    delete(hSlice);title('');delete(hc);
end
%break;
end


% --------------------------------------------------------------------
function varargout = pushbutton4_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton4.
disp('pushbutton4 Callback not implemented yet.')

set(h,'userdata',1);



% --------------------------------------------------------------------
function varargout = popupmenu1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.popupmenu1.
disp('popupmenu1 Callback not implemented yet.')



% --------------------------------------------------------------------
function varargout = slider1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.slider1.
disp('slider1 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = slider2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.slider2.
disp('slider2 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = slider3_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.slider3.
disp('slider3 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = slider4_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.slider4.
disp('slider4 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = slider5_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.slider5.
disp('slider5 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = popupmenu2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.popupmenu2.
disp('popupmenu2 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = pushbutton5_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton5.
disp('pushbutton5 Callback not implemented yet.')

% Create a slice
UD_push2=get(handles.pushbutton2,'userdata');
if(~length(UD_push2))
    uiwait(errordlg('No data','error','modal'));
end

longSlice = UD_push2(1).slongr;
latSlice  = [];
xSlice    = [];
str11=get(handles.popupmenu2,'string');
% Make a new one
if(isstr(str11))
    str1{1} = str11;
else
    str1 = str11;
end
str1{length(str1)+1} = ['slice',num2str(length(str1))];
set(handles.popupmenu2,'string',str1);

% --------------------------------------------------------------------
function varargout = pushbutton6_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton6.
disp('pushbutton6 Callback not implemented yet.')

str1=get(handles.popupmenu2,'string');
% Remove one
if(length(str1)==1 | isstr(str1))
    uiwait(errordlg('No slices in memory','error','modal'));
    return;
end
str1 = str1(1:end-1);
set(handles.popupmenu2,'string',str1);
