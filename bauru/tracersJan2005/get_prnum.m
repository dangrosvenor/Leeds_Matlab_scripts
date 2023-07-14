function imphys=get_prnum(dgcell)

dgs={'PGDEP','PGMLT', ...  %1-2
          'PRAUT',   'PGSHD',   'PRACW',   'PSMLT', ... %3-6
          'PIMLT',   'PSAUT',   'PSDEP',   'PIACR_G', ... %7-10
          'PSACI',   'PGACI',   'PGACW',   'PGACS', ... %11-14
          'PGACR',   'PSACR',   'PRACS',   'PSACW', ... %15-18
          'PGAUT',   'PGFR',    'PRACI_G', 'PGWET', ... %19-22
          'PGDRY',   'PGSUB',   'PSSUB',   'PREVP', ... %23-26
          'PISUB',   'DQI  ',   'PIHAL',   'PIPRM', ... %27-30
          'PIDEP',   'PIACW',   'PICNT',   'PIFRW', ... %31-34
          'PIACR_S', 'PRACI_S', 'PRACI',   'PIACR', ... %35-38
          'RIACR_S', 'RSACR',   'RIACR_G', 'RGFR', ...  %39-42
          'RGACS',   'RGACR',   'RSAUT',   'RIACI', ... %43-46
          'RSACS',   'RSBRK'  ...                       %47-48
                };
            
for idg=1:length(dgcell)            
	for imp=1:length(dgs)
		if strcmp(dgs{imp},dgcell{idg})==1
            imphys(idg)=imp;
            break
		end
	end
end