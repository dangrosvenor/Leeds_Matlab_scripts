function hh = wind_barbs(varargin)
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
%   H = QUIVER(...) returns a vector of line handles.
%
%   Example:
%      [x,y] = meshgrid(-2:.2:2,-1:.15:1);
%      z = x .* exp(-x.^2 - y.^2); [px,py] = gradient(z,.2,.15);
%      contour(x,y,z), hold on
%      quiver(x,y,px,py), hold off, axis image
%
%   See also FEATHER, QUIVER3, PLOT.

%   Clay M. Thompson 3-3-94
%   Copyright 1984-2000 The MathWorks, Inc. 
%   $Revision: 5.19 $  $Date: 2000/06/02 04:30:52 $



% Arrow head parameters
alpha = 0.1; % Size of arrow head relative to the length of the vector
beta = 3;  % Width of the base of the arrow head relative to the length
autoscale = 1; % Autoscale if ~= 0 then scale by this.
plotarrows = 1; % Plot arrows
sym = '';

filled = 0;
ls = '-';
ms = '';
col = '';

nin = nargin;
% Parse the string inputs
while isstr(varargin{nin}),
  vv = varargin{nin};
  if ~isempty(vv) & strcmp(lower(vv(1)),'f')
    filled = 1;
    nin = nin-1;
  else
    [l,c,m,msg] = colstyle(vv);
    if ~isempty(msg), 
      error(sprintf('Unknown option "%s".',vv));
    end
    if ~isempty(l), ls = l; end
    if ~isempty(c), col = c; end
    if ~isempty(m), ms = m; plotarrows = 1; end
    if isequal(m,'.'), ms = ''; end % Don't plot '.'
    nin = nin-1;
  end
end

error(nargchk(2,5,nin));

% Check numeric input arguments
if nin<4, % quiver(u,v) or quiver(u,v,s)
  [msg,x,y,u,v] = xyzchk(varargin{1:2});
else
  [msg,x,y,u,v] = xyzchk(varargin{1:4});
end
if ~isempty(msg), error(msg); end

if nin==3 | nin==5, % quiver(u,v,s) or quiver(x,y,u,v,s)
  autoscale = varargin{nin};
end
% Scalar expand u,v
if prod(size(u))==1, u = u(ones(size(x))); end
if prod(size(v))==1, v = v(ones(size(u))); end
% change to direction vectors
uu_orig=u;vv_orig=v;
u=u./sqrt(u.^2+v.^2);
v=v./sqrt(u.^2+v.^2);

if autoscale,
  % Base autoscale value on average spacing in the x and y
  % directions.  Estimate number of points in each direction as
  % either the size of the input arrays or the effective square
  % spacing if x and y are vectors.
  if min(size(x))==1, n=sqrt(prod(size(x))); m=n; else [m,n]=size(x); end
  %delx = diff([min(x(:)) max(x(:))])/n;
  %dely = diff([min(y(:)) max(y(:))])/m;
  delx = min(diff(x(:)));
  dely = min(diff(y(:)));
  del = delx.^2 + dely.^2;

  if del>0
    len = sqrt((u.^2 + v.^2)/del);
    maxlen = max(len(:));
  else
    maxlen = 0;
  end
  
  if maxlen>0
    autoscale = autoscale*0.9 / maxlen;
  else
    autoscale = autoscale*0.9;
  end
  u = u*autoscale; v = v*autoscale;
end

ax = newplot;
next = lower(get(ax,'NextPlot'));
hold_state = ishold;

% Make velocity vectors
x = x(:).'; y = y(:).';
u = u(:).'; v = v(:).';
uu_orig = uu_orig(:).'; vv_orig = vv_orig(:).';
uu = [x;x+u;repmat(NaN,size(u))];% each column is a vector
vv = [y;y+v;repmat(NaN,size(u))];
h1 = plot(uu(:),vv(:),'k','linewidth',2);

% this is where the wind barbs go
if plotarrows,
  % Make arrow heads and plot them
  % treat the 50m/s barbs as separate
  
  % treat 50m/s, 10m/s and 5m/s
  % 50m/s
  ind_50s=zeros(size(uu_orig));
  index_50s = find(sqrt(uu_orig.^2+vv_orig.^2)>=50);
  ind_50s(index_50s)=floor(sqrt(uu_orig(index_50s).^2+vv_orig(index_50s).^2)./50)
  % 10m/s
  ind_10s=zeros(size(uu_orig));
  index_10s = find((sqrt(uu_orig.^2+vv_orig.^2)-ind_50s.*50)>=10);
  ind_10s(index_10s)=...
      floor((sqrt(uu_orig(index_10s).^2+vv_orig(index_10s).^2)-ind_50s(index_10s).*50)./10);
  % 5m/s
  ind_5s=zeros(size(uu_orig));
  % everything <10
  index_5s = find((sqrt(uu_orig.^2+vv_orig.^2)-ind_50s.*50-ind_10s.*10)>=0);
  ind_5s(index_5s)= ...
      floor((sqrt(uu_orig(index_5s).^2+vv_orig(index_5s).^2)-ind_50s(index_5s).*50-ind_10s(index_5s).*10)./5);
  % number of loops etc
  max_50=max(ind_50s);
  %ind_50s(:)=0;
  max_10=max(ind_10s);
  max_5=max(ind_5s);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %data=local_tr_ep(x_par,w_par,w_perp,num,len_fac,flag)
  %data=local_tr_sp(x_par,w_par,num,flag)
  %flag_ind_u=zeros(size(u));
  %flag_ind_v=zeros(size(v));
  
  % treating 50s separate - different separator
  for ic=1:max_50
      %flag_ind_u(:)=0;
      %flag_ind_v(:)=0;
      flag=find(ind_50s<ic);
      hu50((1+3*(ic-1)):(3+3*(ic-1)),:)=...
          [local_tr_sp(x,u,ic-1,flag);...
              local_tr_ep(x,u,v,ic-1,1,flag,alpha,beta);repmat(NaN,size(u))];
      hv50((1+3*(ic-1)):(3+3*(ic-1)),:)=...
          [local_tr_sp(y,v,ic-1,flag); ...
              local_tr_ep(y,v,u,ic-1,1,flag,alpha,-beta);repmat(NaN,size(v))];
  end
  % treating 10s separate - different separator
  for ic=1:max_10
      %flag_ind_u(:)=0;
      %flag_ind_v(:)=0;
      flag=find(ind_10s<ic);
      hu10((1+3*(ic-1)):(3+3*(ic-1)),:)=...
          [local_tr_sp(x,u,ind_50s+ic-1,flag);...
              local_tr_ep(x,u,v,ind_50s+ic-1,1,flag,alpha,beta);repmat(NaN,size(u))];% problem
      hv10((1+3*(ic-1)):(3+3*(ic-1)),:)=...
          [local_tr_sp(y,v,ind_50s+ic-1,flag); ...
              local_tr_ep(y,v,u,ind_50s+ic-1,1,flag,alpha,-beta);repmat(NaN,size(v))];
  end
  % treating 5s separate - different separator
  for ic=1:max_5
      %flag_ind_u(:)=0;
      %flag_ind_v(:)=0;
      flag=find(ind_5s<ic);
      hu5((1+3*(ic-1)):(3+3*(ic-1)),:)=...
          [local_tr_sp(x,u,ind_50s+ind_10s+ic-1,flag);...
              local_tr_ep(x,u,v,ind_50s+ind_10s+ic-1,0.5,flag,alpha,beta);repmat(NaN,size(u))];% problem
      hv5((1+3*(ic-1)):(3+3*(ic-1)),:)=...
          [local_tr_sp(y,v,ind_50s+ind_10s+ic-1,flag); ...
              local_tr_ep(y,v,u,ind_50s+ind_10s+ic-1,0.5,flag,alpha,-beta);repmat(NaN,size(v))];
  end
  %hu = [x+u;x+u+alpha*(u-beta*(v+eps));repmat(NaN,size(u)); ...
%x+u-0.1*u;x+u-0.1*u+alpha*(u-0.1*u-beta*(v+eps));repmat(NaN,size(u)); ...
%x+u-0.2*u;x+u-0.2*u+0.5.*alpha*(u-0.2*u-beta*(v+eps));repmat(NaN,size(u))];
 % hv = [y+v;y+v+alpha*(v+beta*(u+eps));repmat(NaN,size(v)); ...
%y+v-0.1*v;y+v-0.1*v+alpha*(v-0.1*v+beta*(u+eps));repmat(NaN,size(v)); ...
%y+v-0.2*v;y+v-0.2*v+0.5.*alpha*(v-0.2*v+beta*(u+eps));repmat(NaN,size(v))];


  %hu = [x+u-alpha*(u+beta*(v+eps));x+u; ...
  %      x+u-alpha*(u-beta*(v+eps));repmat(NaN,size(u))];
  %hv = [y+v-alpha*(v-beta*(u+eps));y+v; ...
  %      y+v-alpha*(v+beta*(u+eps));repmat(NaN,size(v))];
  hold on
  % now do plots %%%%%%%%%%%%%%%%
  if(exist('hu50','var'))
    h2 = plot(hu50(:),hv50(:),'r','linewidth',2);
  else
    h2=[];
  end
  if(exist('hu10','var'))
    %h4 = plot(hu10(:),hv10(:),[col ls]);
    h4 = plot(hu10(:),hv10(:),'k','linewidth',2);
  else
    h4=[];
  end
  if(exist('hu5','var'))
    h5 = plot(hu5(:),hv5(:),'k','linewidth',2);
  else
    h5=[];
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
  h2 = [];
  h4 = [];
  h5 = [];
end

if ~isempty(ms), % Plot marker on base
  hu = x; hv = y;
  hold on
  h3 = plot(hu(:),hv(:),'ok','linewidth',2);
  if filled, set(h3,'markerfacecolor',get(h1,'color')); end
else
  h3 = [];
end

if ~hold_state, hold off, view(2); set(ax,'NextPlot',next); end

if nargout>0, hh = [h1;h2;h3]; end

function data=local_tr_ep(x_par,w_par,w_perp,num,len_fac,flag,alpha,beta)
% tranform coordinates for end point

%global distance;
data=x_par+w_par-0.2.*num.*w_par+len_fac.*alpha* ...
    (w_par-0.2.*num.*w_par-beta.*(w_perp+eps));
data(flag)=nan;
%x+u-0.2*u;x+u-0.2*u+0.5.*alpha*(u-0.2*u-beta*(v+eps));repmat(NaN,size(u))];

function data=local_tr_sp(x_par,w_par,num,flag)
% tranform coordinates for start point
global alpha;
global beta;
%global distance;
data=x_par+w_par-0.2.*num.*w_par;
data(flag)=nan;
