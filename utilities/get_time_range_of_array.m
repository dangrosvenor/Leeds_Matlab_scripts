function [out, time_out, time_inds, dtime_match] = get_time_range_of_array(array_in,time_matlab,time_choice,dim)
% function [out, time_out, time_inds] = get_time_range_of_array(array_in,time_matlab,time_choice,dim)
% Finds the correct portion of an array given a set of time ranges or
% specific times and some tolerances
% To use the time range method specify time_matlab.time_range to be an [M 2] 
% sized array where have M different ranges within time_matlab
% It will return all values between those ranges.
% For specific times need to remove time_choice.time_range and then specify
% time_choice.time_specific and time_choice.tol for the times required and
% the tolerances.
% If just want the time indies then set array_in=[]
% If just want the closest match in time then set time_choice.find_nearest=1
% Could just use interp1 with 'nearest' option?

if isfield(time_choice,'find_nearest')==0
    time_choice.find_nearest=0;
end

if ndims(array_in)<3
    dtime_match=0;
else
    dtime_match=zeros(size(array_in));
end
    


%Sort out whether requesting a range or specific times
if isfield(time_choice,'time_range')==1
    time_range = time_choice.time_range;
    if size(time_range,2) > 2
        error('Time range should be of size [M 2]');
    end
    time_inds=[];
    for irange=1:size(time_range,1)
        time_inds2 = find(time_matlab >= time_choice.time_range(irange,1) & time_matlab <= time_choice.time_range(irange,2) );
        time_inds = cat(2,time_inds,time_inds2);
    end
%  time_inds=[1 3 5];


else
    time_specific = time_choice.time_specific;
    nT = length(time_specific);
    if ndims(array_in)<3
        dtime_match=0;
    else
        dtime_match=zeros([size(array_in,1) size(array_in,2) nT]); %fix to work with dim
    end

    for it=1:nT
        time = time_specific(it);
        dtime = abs(time_matlab - time);
        if time_choice.find_nearest==0
            
            if ndims(time_matlab)<3
                
                ifind = find(  dtime < time_choice.tol );
                
                if length(ifind)==0
                    error('Could not find the specified times within the chosen tolerance.');
                elseif length(ifind)>1
                    error('More than one time found within the chosen tolerance.');
                else
                    time_inds(it)=ifind;
                end
                
                
            else
                
                %[dtime_match(:,:,it),time_inds(:,:,it)]=min(dtime,[],dim);
                %so here time_inds are the time index for the matching time
                %for every location in the 2D lat/long array
                
                [time_inds(:,:,it),out(:,:,it),min_vals] = min_column_inds(permute(dtime,[3 1 2]), permute(array_in,[3 1 2]));
                %min_column_inds returns the values of the second array (arr2) at the min of the first
                % (arr) over the 1st dimension. Also returns the indices (inds).
                [time_out(:,:,it)] = apply_inds_from_one_dim_in_3d_array(time_inds(:,:,it), permute(time_matlab,[3 1 2]));
                %[time_inds(:,:,it),time_out(:,:,it)] = min_column_inds(permute(dtime,[3 1 2]), permute(time_matlab,[3 1 2]));
                
                inan=find(abs(min_vals)>time_choice.tol);
                temp=out(:,:,it);
                temp(inan)=NaN;
                out(:,:,it)=temp;
                
                
            end

        else
            %Otherwise find the closest matching time
            if ndims(time_matlab)<3
                [dtime_match(it),time_inds(it)]=min(dtime);               
            end
                       
        end
             
    end   
    
    
end


if ndims(time_matlab)<3
    
    time_out = time_matlab(time_inds(:));    
    
    if size(array_in,1)~=0
        
        %Make dim the first dimension
        siz=size(array_in); siz2=siz;
        inds = 1:length(siz); inds2=inds;
        inds2(dim)=[];
        siz2(dim)=[];
        per_inds = [dim inds2];
        array_in = permute(array_in,per_inds);
        
        %Figure out how to get the original array back after permuting...
        %  Actaully can just use ipermute (inverse permute) with orig indices.
        % per_inds_out = NaN*ones(size(inds));
        % a = 2:length(siz);
        % i = inds;
        % i(dim)=[];
        % per_inds_out(inds2) = a;
        % per_inds_out(dim) = 1;                
        
        out = array_in(time_inds,:);
        out=reshape(out,[length(time_inds) siz2]);
                       
        
        % Permute back to original order
        out  = ipermute(out,[per_inds]);
    else
        out=NaN;
    end
    
end

