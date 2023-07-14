function i=integrate(x,v)

i=0.5*sum( ( v(1:end-1)+v(2:end) ) .* ( x(2:end)-x(1:end-1) ) );