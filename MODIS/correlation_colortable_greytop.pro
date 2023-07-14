pro correlation_colortable_greytop,perc,greylev

;perc=fraction of full scale in centre that should be grey 
;greylev=grey brightness (255 white, 0 black)

loadct,33
tvlct,r,g,b,/Get

red=r
grn=g
blu=b

;squash 0:154 into 0:128-perc*0.5*256
nx=128-fix(perc*0.5*256)
xint=135.0*findgen(128-fix(perc*0.5*256))/(128.-perc*0.5*256)
x=findgen(256)

redneg=interpol(float(red),x,xint)
grnneg=interpol(float(grn),x,xint)
bluneg=interpol(float(blu),x,xint)

xint=154.+101.0*findgen(128-fix(perc*0.5*256))/(128.-perc*0.5*256)

redpos=interpol(float(red),x,xint)
grnpos=interpol(float(grn),x,xint)
blupos=interpol(float(blu),x,xint)

nsofar=256-n_elements([redneg,redpos])
if nsofar gt 0 then begin
  mid=bytarr(nsofar)
  red=byte([redneg,mid,redpos])
  grn=byte([grnneg,mid,grnpos])
  blu=byte([bluneg,mid,blupos])
;  red(128-fix(perc*256.0*0.5):128+fix(perc*256*0.5))=byte(greylev)
;  grn(128-fix(perc*256.0*0.5):128+fix(perc*256*0.5))=byte(greylev)
;  blu(128-fix(perc*256.0*0.5):128+fix(perc*256*0.5))=byte(greylev)
endif else begin
  red=byte([redneg,redpos])
  grn=byte([grnneg,grnpos])
  blu=byte([bluneg,blupos])   
endelse


;assign 255 to white
red(255)=255
grn(255)=255
blu(255)=255

;assign 0 to black
red(0)=000
grn(0)=000
blu(0)=000

;assign 254 to grey
red(254)=byte(greylev)
grn(254)=byte(greylev)
blu(254)=byte(greylev)

tvlct,red,grn,blu



end
