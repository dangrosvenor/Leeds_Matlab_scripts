function rel=calc_relAz(sens,sol)
%calculate relative azimuth angle from sensor and solar Az

%Relative azimuth is defined as zero when the sensor is looking
            %into the Sun - as if there as no scattering (just
%             %transmission). Backscatter is then at 180 degrees. So do the
%             difference of sensor and solar and subtract 180 so that if
%             the difference is 180 (forward scatter) then will get relAZ=0
%NOTE - %with Joint files only get the relative AZ, but have kept the variable
%Solar_Azimuth_Angle for consistency with full L2 and have set the sensorAZ to be the
%relative AZ. But are the MODIS Relative_Azimuth_Angles with the 180
%already removed? I.e. do they follow the above convention? Check using L2?

rel = abs(sens - sol - 180);
rel(rel>180) = abs( 360-rel(rel>180) );

