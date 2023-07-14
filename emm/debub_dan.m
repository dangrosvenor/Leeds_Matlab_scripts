exname='g:/emm_runs/lwc_thresh2/debug_dan';

debug=dlmread(exname,' ');

for j=1:3
	i=find(debug(:,1)==j);
	debug2(j).d=debug(i,:);
end