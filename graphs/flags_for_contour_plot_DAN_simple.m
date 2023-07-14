clear h min

iabc=1; %flag to tell it to add (a), (b), (c), etc. to runName(i).nam
comp='uni';
iovride_conts=0;

isquare=0; %flag to make axes square

iadd_overlay = 0;
icont_neg=0;
idraw_streamlines=0;
right_side_extra_bits=1; %switch to put the extreme values in blocks on right hand side of the plot
dan_test=0; %test bit of script, 28th Oct 2009

if ~exist('iplotselect'); iplotselect=0; end
if ~exist('bigcbar'); bigcbar=0; end
if ~exist('lememm'); lememm=0; end
if ~exist('iabc'); iabc=0; end
if ~exist('comp'); comp='uni'; end
if ~exist('icont_extra'); icont_extra=0; end
if ~exist('i3d'); i3d=0; end
if ~exist('wrap2d'); wrap2d=0; end

if iplotselect==1
    plotcase=plotcases(iplot);
end

clear zz timesTH diff

if ~exist('idir'); idir=1; end
if ~exist('ieps'); ieps=0; end
if ~exist('isamescale'); isamescale=0; end
if ~exist('subplotting'); subplotting=0; end
if ~exist('onexlabel'); onexlabel=0; end      %flag to make it so that only the bottom plot has the xlabel on it
if ~exist('iplot'); iplot=1; end

if exist('npess2')
    if plotcase~=65
        npes=npess2(idir);
    end
end

%plots time height plots

ilem=1;
icont=1;

if ieps==1
    fsize=12;
elseif subplotting==1
    fsize=18;
else
    fsize=14; %26 best for paper figures when more than one plot is being used
end

isave=0;
%load c:/matlabr12/work/bauru/casestudy/forcecons/diag/profnew+ccn_2-169;

icolmap=1; %flag to set to another colormap (defined by cmap)

iutc=1; %flag for time axis to be labelled as UTC (otherwise is labelled Local Time)
%add_ground_height=0.62; %height to add to the vertical axis to account for level of ground abv msl.
%add_ground_height=0; %height to add to the vertical axis to account for level of ground abv msl.
%add_ground_height=1.0; %height to add to the vertical axis to account for level of ground abv msl.


minZ=0e3;
maxZ=25e3;  %19000;
maxtr=1.0;
%timesTH=Time;

vectorf=0;

jmax=5; %max no. plots on one screen

a1=1;
a2=2; %values for subplot(ai,b,i)

izmin=2;
izmax=2;
f=1e6*28.97/18;

%iminovr=zeros([1 10]);
%imaxovr=zeros([1 10]);
notsame=0; %flag to plot each plot individually in terms of colour scale
offset=0;
clines=0; % when zero removes black contour lines
nocbar=0;
sig=2; %no. sig figs for contour values
nplots2d=1;
itimestamp=0;
manclab=0; %flag to set contour labels manually with mouse
idirstamp=0; %flag to write the directory of the results in the corner

normcbar=0;

dumpint=300; %dump interval in seconds

i2d=0;
izovr=0; %flag to say are setting own z axis
itimelab=0; %flag to say that x axis is time and so should be restricted to <24
figlab='2d contour plot';
iylim=0; %flag to override the setting of the y part of axis (e.g. so can have axis larger than that plotted)
ixlim=0;

nplots2d2d=1;

ncont=25;
clab=1; %flag to label contours
i=1;

phase=1;


clear max timestxt pdat minC maxC

logflag=0;
dlogflag=0;




%%%%%%%%%%%%%%%%%%%%%%%%    start of old plot case 577 %%%%%%%%%%%%%%%%%%%%%
idirstamp=1;


fact=1e6*28.97/18;
logflag=0;
itimestamp=1;




imaxovr=0;
iminovr=0;

clines=1; %makes black contour lines appear
clab=0;

ilem=0;


i2d=2; %tells it to label x axis in km



%%%%%%%%%%%%%%%%%%%%%%%%    case 'wrf_plot'

