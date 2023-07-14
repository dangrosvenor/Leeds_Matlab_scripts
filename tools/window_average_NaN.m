function [x,y,N_notnan] = window_average_NaN(xdat,ydat,Nwindow,operation_flag,nthresh)
%function [x,y] = window_average(xdat,ydat,Nwindow,operation_flag)
%For an Nd array it assumes that the FIRST dimension is the time dimension
%- others are collapsed and restored after.

if nargin>4
   i_use_Nthresh=1;
else
   i_use_Nthresh=0;
end

siz=size(ydat);
nd=length(find(siz>1)); %Number of non-singleton dimensions
if nd>1
    nT=size(ydat,1);
    ydat=ydat(:,:); %collapse into 2 dim array
else
    nT=length(ydat);    
end

ntot = nT-Nwindow+1;

%Pre-alocate arrays for increased speed
y=NaN*ones([ntot size(ydat,2)]);
N_notnan=NaN*ones([ntot size(ydat,2)]);


ntenth = round(ntot/10);
fprintf(1,'\n%d: ', ntot);
icount=0;
for istart=1:ntot
    icount=icount+1;
    if icount==ntenth
        fprintf(1,'%d ',istart);
        icount=0;
    end
    
    if nd>1
        [d,n] = meanNoNan(ydat(istart:istart+Nwindow-1,:),1,operation_flag);   
        if i_use_Nthresh==1 %NaN data where fewer than nthresh non_NaN points go into the mean
            inan = find(n < nthresh);
            d(inan) = NaN;
        end
%         if maxALL(n)==343
%             'yes'
%         end
        y(istart,:) = d;        
        N_notnan(istart,:) = n;
    else
        [d,n] = meanNoNan(ydat(istart:istart+Nwindow-1),1,operation_flag);
        if i_use_Nthresh==1 %NaN data where fewer than nthresh non_NaN points go into the mean
            inan = find(n < nthresh);
            d(inan) = NaN;
        end
        y(istart) = d;
        N_notnan(istart) = n;
    end    
end

if nd>1
    y=reshape(y,[size(y,1) siz(2:end)]);
    N_notnan=reshape(N_notnan,[size(y,1) siz(2:end)]);    
end
    



%now fix x so that it corresponds to the mid-points of the data over which
%the window was applied
if mod(Nwindow,2)==0  %even number - in this case the window mid-points falls between datapoints in x
    %and so we will find all of the points halfway between 
    x =0.5*(xdat(Nwindow/2:end-Nwindow/2) + xdat(Nwindow/2+1:end-Nwindow/2+1));
else %for odd Nwindow values we can just use the exact x datapoints starting at the mid-point of the first window
    x = xdat(Nwindow/2+0.5:end-Nwindow/2+0.5);
end

