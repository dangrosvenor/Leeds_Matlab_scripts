%Carry through weights from last week
weights(:,iweek) = weights(:,iweek-1);

dat_week = choices(:,:,iweek); %5*N_pers
%Find majority choice after weighting
weights_week = ( repmat(weights(:,iweek),[1 5]) )';
% But if not using the weights as done previously then set to one.
if maj_system==1
    weights_week = ones(size(weights_week));
end

tot = meanNoNan(dat_week.*weights_week,2,'sum');
tot_unweighted = meanNoNan(dat_week,2,'sum');
[maj_val,imaj] = max(tot); %imaj becomes the index of the majority day.


fprintf(1,'\n---- Week %d -----\n',iweek);
fprintf(1,'Day           Total votes     Weighted score\n');
for i=1:length(days)
    fprintf(1,'%s      %f        %f \n',days{i},tot_unweighted(i),tot(i));
end


%Check whether two days got the same score
irep = find(tot==maj_val);
if length(irep)>1
    if choice_override==0  %The unset value
        error('Scores indicate a tie - need to enter which day was chosen!');
    else
        imaj=choice_override;
        choice_override=0;
    end
end

fprintf(1,'\nDay chosen = %s\n',days{imaj});
fprintf(1,'\n');

%Now find people who could not do the chosen day
ipers_denied = find(dat_week(imaj,:)==0);

%Increase the weights by 0.25 of those people
weights(ipers_denied,iweek) = weights(ipers_denied,iweek)+0.25;

%Decrease the weights by 0.25 of those people who got their choice
ipers = find(dat_week(imaj,:)==1);
weights(ipers,iweek) = weights(ipers,iweek)-0.25;
%But limit to >1
weights(weights<1)=1;



fprintf(1,'Name     Old weights     New weights\n');
for i=1:length(peeps)
    fprintf(1,'%s      %f      %f\n',peeps{i},weights_week(1,i), weights(i,iweek));
end
fprintf(1,'\n');

