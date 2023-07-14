function v=vari(mean,vals,nps)
	tot=sum(nps);
	sqdiff=(mean-vals).^2;
	v=sum(nps.*sqdiff(1:length(nps)))./(tot-1);
