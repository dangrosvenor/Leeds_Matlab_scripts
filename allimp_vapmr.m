if strmatch('Vap MR & ice NC',constr)
            f=1e6*28.97/18;
            
            if (usediags==1)
                timei=mod(GridDan(idir).t(jj)+3,24);
            else
                timei=mod(SER(end,1)/3600 + 19.67,24);
            end
            
            timlab=num2str(round2(timei,2));
            noplot=1;
            output=0;
            
            iplot57=1;
            for i57=[8]
                switch i57
                    case 1
                        %exdirA=[direcDan(jc).dir 'results/totMR+iceNC/'];
                        exdirA=[direcDan(jc).dir 'results/totMRzoom/'];
                        iplot57=1;
                        i577='totwater';
                    case 2
                        exdirA=[direcDan(jc).dir 'results/vapMR+iceNC/'];
                    case 3
                        i577='tot_condensate';
                        exdirA=[direcDan(jc).dir 'results/totcond/'];
                        iplot57=1;
                        %exdirA=['c:/documents and settings/G/my documents/lem/totCondzoom/'];
                    case 4
                        i577='si';   
                        exdirA=[direcDan(jc).dir 'results/icesupersatPlume/'];
                        exdirA=['c:/documents and settings/G/my documents/lem/icesupersatPlume/'];
                        iplot57=1;
                    case 5
                        i577='temppert';                
                        exdirA=[direcDan(jc).dir 'results/temp_pert_emm_thermal/']; 
                    case 6
                        i577='rhopert';                
                        exdirA=[direcDan(jc).dir 'results/rhopert/']; 
                        iplot57=0; %just want to calculate data without plotting
                    case 7
                        i577='hydbal';                
                        exdirA=[direcDan(jc).dir 'results/hydbal/'];
                    case 8
                        exdirA=[direcDan(jc).dir 'results/vapMR/'];
                   %     exdirA=['c:/documents and settings/G/my documents/lem/1km_res/iceMR_250m1000km/'];

                        iplot57=1;
                        i577='vapour';  
                    case 9
                        exdirA=[direcDan(jc).dir 'results/LNB/'];
                        iplot57=0;
                        i577='lnb';   
                    case 10
                        exdirA=[direcDan(jc).dir 'results/fall/'];
                        iplot57=0;
                        i577='fall';   
                    case 11
                        exdirA=[direcDan(jc).dir 'results/potemp/'];
                       % exdirA=['c:/documents and settings/G/my documents/lem/1km_res/potemp2/'];
                        i577='potemp';  
                     case 12
                        iplot57=0;
                        i577='temps';   
                      case 13
                        exdirA=[direcDan(jc).dir 'results/vertvel_1-250km/'];
                        iplot57=1;
                        i577='vertvel';  
                      case 14
                        iplot57=0;
                        i577='vapdists';    
                    case 15
                       i577='icesatMR';    
                        exdirA=[direcDan(jc).dir 'results/icesatMR_1-250km/'];
                        iplot57=1;
                    case 16
                        iplot57=0;
                        i577='temppert_signchange_timser_LNB';
                    case 17
                        iplot57=0;
                        i577='lowtracer';
                        exdirA=[direcDan(jc).dir 'results/low_tracer/'];
                    case 18
                        iplot57=1;
                        i577='inc';
                        exdirA=[direcDan(jc).dir 'results/inc/'];
                    case 19
                        iplot57=0;
                        i577='meaninc';           
                        
                    end
                
                if iplot57==1
                    bigcbar=0; %flag for one big colorbar instead of one for each subplot
					isamescale2=0;
					subplotting=0; %flag for whether to use subplots instead of new figures
                    lememm=0;
                    
					plotTimeHeightVap3;
                    jc=1;
					%exname=strcat(exdirA,int2str(jj),'-',int2str(k),'-',textdataDan(jc).text,'-',int2str(cono),'.emf')
                    exname=strcat(exdirA,int2str(jj),'-',int2str(3),'-',textdataDan(jc).text,'-',int2str(cono),'.emf');
					set(hf,'paperpositionmode','auto');
					%print(gcf,'-djpeg','-r350',exname);
					print(hf,'-dmeta',exname);
					close(hf);
                else
                    calcsfor2dplots
                end
                    
                    
                    
            end
            
          
     end