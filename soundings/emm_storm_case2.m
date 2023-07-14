ps_Pa(1) =  595.84 ; temp_K(1) =   0.34; td_K(1) = -0.95;		zRead(1) = 5.;
ps_Pa(2 ) =  553.60; temp_K(2) =  -2.48; td_K(2) =  -11.14;		zRead(2) = 6.;
ps_Pa( 3) =  542.76; temp_K(3) =  -3.83; td_K(3) = -12.76;		zRead(3) = 6.2;
ps_Pa( 4) =  505.85; temp_K(4) =   -6.31; td_K(4) = -15.47;		zRead(4) = 6.6;
ps_Pa( 5) = 481.29; temp_K(5) =  -10.00; td_K(5) =  -14.39;		zRead(5) = 6.8;
ps_Pa( 6) =  450.45; temp_K(6) =  -12.60; td_K(6) =  -18.57;	zRead(6) = 7.5;
ps_Pa( 7) = 425.30; temp_K(7) =  -14.89; td_K(7) = -22.73;		zRead(7) = 7.8;
ps_Pa( 8) = 406.12  ; temp_K(8) =  -16.91; td_K(8) = -23.54 ;	zRead(8) = 8.;
ps_Pa( 9) =  375.00; temp_K(9) =  -21.10; td_K(9) = -27.10;		zRead(9) = 11.;
ps_Pa( 10) =  300.00; temp_K(10) =  -32.90; td_K(10) =  -62.90;	zRead(10) = 12.;
ps_Pa( 11) =  267.00; temp_K(11) =  -40.10; td_K(11) =  -52.10;	zRead(11) = 14.;

for i=12:-1:2
    ps_Pa(i)=ps_Pa(i-1);
    temp_K(i)=temp_K(i-1);
    td_K(i)=td_K(i-1);
    zRead(i)=zRead(i-1);
end

ps_Pa(1) = 642.86; temp_K(1) = 5.77; td_K(1) = 1.44;			zRead(1) = 3.;

    

