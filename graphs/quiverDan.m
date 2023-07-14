function hh = quiver(varargin)
%QUIVER Quiver plot.
%   QUIVER(X,Y,U,V) plots velocity vectors as arrows with components (u,v)
%   at the points (x,y).  The matrices X,Y,U,V must all be the same size
%   and contain corresponding position and velocity components (X and Y
%   can also be vectors to specify a uniform grid).  QUIVER automatically
%   scales the arrows to fit within the grid.
%
%   QUIVER(U,V) plots velocity vectors at equally spaced points in
%   the x-y plane.
%
%   QUIVER(U,V,S) or QUIVER(X,Y,U,V,S) automatically scales the 
%   arrows to fit within the grid and then stretches them by S.  Use
%   S=0 to plot the arrows without the automatic scaling.
%
%   QUIVER(...,LINESPEC) uses the plot linestyle specified for
%   the velocity vectors.  Any marker in LINESPEC is drawn at the base
%   instead of an arrow on the tip.  Use a marker of '.' to specify
%   no marker at all.  See PLOT for other possibilities.
%
%   QUIVER(...,'filled') fills any markers specified.
%
%   QUIVER(AX,...) plots into AX instead of GCA.
%
%   H = QUIVER(...) returns a quivergroup handle.
%
%   Example:
%      [x,y] = meshgrid(-2:.2:2,-1:.15:1);
%      z = x .* exp(-x.^2 - y.^2); [px,py] = gradient(z,.2,.15);
%      contour(x,y,z), hold on
%      quiver(x,y,px,py), hold off, axis image
%
%   See also FEATHER, QUIVER3, PLOT.

%   Backwards compatibility
%   QUIVER('v6',...) creates line objects instead of a quivergroup
%   object for compatibility with MATLAB 6.5 and earlier.

%   Clay M. Thompson 3-3-94
%   Copyright 1984-2006 The MathWorks, Inc. 
%   $Revision: 5.21.6.14 $  $Date: 2006/12/20 07:18:11 $

[v6,args] = usev6plotapi(varargin{:});
if v6
  h = Lquiverv6(args{:});
else
  error(nargchk(1,inf,nargin,'struct'));
  % Parse possible Axes input
  [cax,args] = axescheck(args{:});
  [pvpairs,unused,unused,msg] = parseargs(args);
  error(msg);
  
  if isempty(cax) || isa(handle(cax),'hg.axes')
      cax = newplot(cax);
      parax = cax;
      hold_state = ishold(cax);
  else
      parax = cax;
      cax = ancestor(cax,'Axes');
      hold_state = true;
  end
  [ls,c] = nextstyle(cax);

  h = specgraph.quivergroup('Color',c,'LineStyle',ls,...
                            'parent',parax,pvpairs{:});

  if ~any(strcmpi('color',pvpairs(1:2:end)))
    set(h,'CodeGenColorMode','auto');
  end
  set(h,'refreshmode','auto');
  if ~hold_state, box(cax,'on'); end
  plotdoneevent(cax,h);
  h = double(h);
end

if nargout>0, hh = h; end

function h = Lquiverv6(varargin)
[cax,args,nargs] = axescheck(varargin{:});

% Arrow head parameters
alpha = 0.33; % Size of arrow head relative to the length of the vector
beta = 0.33;  % Width of the base of the arrow head relative to the length
% NOTE - this is changed later on to correct for when the vertical axis is smaller than the horizontal

autoscale = 1; % Autoscale if ~= 0 then scale by this.
plotarrows = 1; % Plot arrows

filled = 0;
ls = '-';
ms = '';
col = '';

nin = nargs;
% Parse the string inputs
while ischar(args{nin}),
  vv = args{nin};
  if ~isempty(vv) && strcmpi(vv(1),'f')
    filled = 1;
    nin = nin-1;
  else
    [l,c,m,msg] = colstyle(vv);
    if ~isempty(msg), 
      error(id('UnknownOption'),'Unknown option "%s".',vv);
    end
    if ~isempty(l), ls = l; end
    if ~isempty(c), col = c; end
    if ~isempty(m), ms = m; plotarrows = 0; end
    if isequal(m,'.'), ms = ''; end % Don't plot '.'
    nin = nin-1;
  end
end

error(nargchk(2,5,nin,'struct'));

% Check numeric input arguments
if nin<4, % quiver(u,v) or quiver(u,v,s)
  [msg,x,y,u,v] = xyzchk(args{1:2});  % "xyzchk(Z,C)"
else
  [msg,x,y,u,v] = xyzchk(args{1:4});  % "xyzchk(X,Y,Z,C)"
end
if ~isempty(msg), 
  if isstruct(msg)
    % make the xyzchk message match quiver's help string:
    msg.message = strrep(msg.message, 'Z', 'U');
    msg.message = strrep(msg.message, 'C', 'V');
    msg.identifier = strrep(msg.identifier, 'Z', 'U');
    msg.identifier = strrep(msg.identifier, 'C', 'V');
  else
    msg = strrep(msg, 'Z', 'U');
    msg = strrep(msg, 'C', 'V');
  end        
  error(msg); %#ok
end

if nin==3 || nin==5, % quiver(u,v,s) or quiver(x,y,u,v,s)
  autoscale = args{nin};
end

% Scalar expand u,v
if numel(u)==1, u = u(ones(size(x))); end
if numel(v)==1, v = v(ones(size(u))); end



%autoscale = 0.05 * autoscale ./ sqrt( (u/(max(x(:))-min(x(:)))/n).^2  + (v/(max(y(:))-min(y(:)))).^2 );
%u=u.*autoscale;
%v=v.*autoscale;


if min(size(x))==1, n=sqrt(numel(x)); m=n; else [m,n]=size(x); end

xmin = min(x(:));
xmax = max(x(:));
ymax = max(y(:));
ymin = min(y(:));

%convert to the co-ordinates ranging from 0 to 1 in x and y (square grid, Dan)
%x = ( x - xmin ) / ( xmax - xmin );
%y = ( y - ymin ) / ( ymax - ymin );


if autoscale,
  % Base autoscale value on average spacing in the x and y
  % directions.  Estimate number of points in each direction as
  % either the size of the input arrays or the effective square
  % spacing if x and y are vectors.
  if min(size(x))==1, n=sqrt(numel(x)); m=n; else [m,n]=size(x); end
  delx = diff([min(x(:)) max(x(:))])/n;
  dely = diff([min(y(:)) max(y(:))])/m;
  del = delx.^2 + dely.^2;
  if del>0
    len = sqrt((u.^2 + v.^2)/del);
    maxlen = max(len(:));
  else
    maxlen = 0;
  end

%autoscale = autoscale * dely/delx;

%autoscale = 0.05 * autoscale ./ sqrt( (u/(max(x(:))-min(x(:)))/n).^2  + (v/(max(y(:))-min(y(:)))).^2 );
  
  if maxlen>0
    autoscale = autoscale*0.9 / maxlen;
  else
    autoscale = autoscale*0.9;
  end
  u = u.*autoscale; v = v.*autoscale;





end

cax = newplot(cax);
next = lower(get(cax,'NextPlot'));
hold_state = ishold(cax);

% Make velocity vectors
x = x(:).'; y = y(:).';
u = u(:).'; v = v(:).';
uu = [x;x+u;repmat(NaN,size(u))];
vv = [y;y+v;repmat(NaN,size(u))];

%uu = xmin + uu*(xmax - xmin);  %to scale back from a square grid to actual x,y coords, Dan
%vv = ymin + vv*(ymax - ymin);

h1 = plot('v6',uu(:),vv(:),[col ls],'parent',cax);

if plotarrows,
  % Make arrow heads and plot them

  delx = diff([min(x(:)) max(x(:))])/n
  dely = diff([min(y(:)) max(y(:))])/m

%beta = beta * 100 * dely/delx %made this change as when the distance covered by the vertical axis was small compared to that of the horizontal the arrow heads were much too large vertically. This scales according to the ratio of the vertical grid size to the horizontal -Dan, 8thOct08

%beta = ones(size(u))*0.33;
%beta = beta .* 
%beta=1000;
  hu = [x+u-alpha*(u+beta*(v+eps));x+u; ...
        x+u-alpha*(u-beta*(v+eps));repmat(NaN,size(u))];
  hv = [y+v-alpha*(v-beta*(u+eps));y+v; ...
        y+v-alpha*(v+beta*(u+eps));repmat(NaN,size(v))];
  hold(cax,'on')

hu = xmin + hu*(xmax - xmin);
hv = ymin + hv*(ymax - ymin);

  h2 = plot('v6',hu(:),hv(:),[col ls],'parent',cax);
else
  h2 = [];
end

if ~isempty(ms), % Plot marker on base
  hu = x; hv = y;
  hold(cax,'on')
  h3 = plot('v6',hu(:),hv(:),[col ms],'parent',cax);
  if filled, set(h3,'markerfacecolor',get(h1,'color')); end
else
  h3 = [];
end

if ~hold_state, hold(cax,'off'), view(cax,2); set(cax,'NextPlot',next); end

h = [h1;h2;h3];

function [pvpairs,args,nargs,msg] = parseargs(args)
% separate pv-pairs from opening arguments
[args,pvpairs] = parseparams(args);
n = 1;
extrapv = {};
% check for 'filled' or LINESPEC
while length(pvpairs) >= 1 && n < 3 && ischar(pvpairs{1})
  arg = lower(pvpairs{1});
  if arg(1) == 'f'
    pvpairs(1) = [];
    extrapv = {'MarkerFaceColor','auto',extrapv{:}};
  else
    [l,c,m,tmsg]=colstyle(pvpairs{1});
    if isempty(tmsg)
      pvpairs(1) = [];
      if ~isempty(l) 
        extrapv = {'linestyle',l,extrapv{:}};
      end
      if ~isempty(c)
        extrapv = {'color',c,extrapv{:}};
      end
      if ~isempty(m)
        extrapv = {'ShowArrowHead','off',extrapv{:}};
        if ~isequal(m,'.')
          extrapv = {'marker',m,extrapv{:}};
        end
      end
    end
  end
  n = n+1;
end
pvpairs = [extrapv pvpairs];
msg = checkpvpairs(pvpairs);
nargs = length(args);
if isa(args{nargs},'double') && (length(args{nargs}) == 1) && ...
      (nargs == 3 || nargs == 5)
  if args{nargs} > 0
    pvpairs = {pvpairs{:},'autoscale','on',...
	       'autoscalefactor',args{nargs}};
  else
    pvpairs = {pvpairs{:},'autoscale','off'};
  end
  nargs = nargs - 1;
end
x = [];
y = [];
u = [];
v = [];
if (nargs == 2)
    u = datachk(args{1});
    v = datachk(args{2});
  pvpairs = {pvpairs{:},'udata',u,'vdata',v};
elseif (nargs == 4)
    x = datachk(args{1});
    y = datachk(args{2});
    u = datachk(args{3});
    v = datachk(args{4});
  pvpairs = {pvpairs{:},'xdata',x,'ydata',y,'udata',u,'vdata',v};
end

if isempty(x)
    % Deal with quiver(U,V) syntax
    if ~isempty(u)
        if ~isequal(size(u),size(v))
            msg.identifier = id('UVSizeMismatch');
            msg.message = 'U and V must be the same size.';
        end
    end
else
    % Deal with quiver(X,Y,U,V) syntax.
    if ~isvector(u)
        msg = xyzcheck(x,y,u,'U');
    else
        % If all four are vectors, xyzcheck is not used as it assumes the
        % third argument should be a matrix.
        if ~isequal(size(u),size(v))
            msg.identifier = id('UVSizeMismatch');
            msg.message = 'U and V must be the same size.';
        elseif ~isequal(size(y),size(u));
            msg.identifier = id('YUSizeMismatch');
            msg.message = 'The size of Y must match the size of U or the number of rows of U.';
        elseif ~isequal(size(x),size(u));
            msg.identifier = id('XUSizeMismatch');
            msg.message = 'The size of X must match the size of U or the number of columns of U.';            
        end
    end
    if isempty(msg) && ~isequal(size(u),size(v))
        msg = [];
        msg.identifier = id('UVSizeMismatch');
        msg.message = 'U and V must be the same size.';
    end
end

function str=id(str)
str = ['MATLAB:quiver:' str];
