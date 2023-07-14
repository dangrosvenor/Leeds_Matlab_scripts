function [time,sza2] = get_time_for_SZA(lat,lon,sza,day)

tstep=4/24;

x0=[day-tstep day+tstep];
root=0;

step=0.2;
Nmax=20;

tstep2=tstep/50;
js=x0(1):tstep2:x0(2);

xa=x0(1);

root=0;
j=0;
while j<length(js) & root==0
    j=j+1;

    xb=js(j);

    sza2=sza;
    i=0;
    while (i<Nmax)

        if sign(sun_pos_fzero(xa,lat,lon,sza2)) == sign(sun_pos_fzero(xb,lat,lon,sza2))
            sza2 = sza2 - step;
            i=i+1;
        else
            root=1;
            i=1e20;
        end

    end

    if root==0

        sza2=sza;
        i=0;
        while (i<Nmax)

            if sign(sun_pos_fzero(xa,lat,lon,sza2)) == sign(sun_pos_fzero(xb,lat,lon,sza2))
                sza2 = sza2 + step;
                i=i+1;
            else
                root=1;
                i=1e20;
            end

        end

    end

end

root
if root==1
    time = fzero(@sun_pos_fzero,[xa xb],[],lat,lon,sza2);
end
