function [fz1,fz2,fz3] = measurements(fx,fxp)
%æ‡¿Î¡ø≤‚
fz1=sqrt((fx(1)-fxp(1))^2+((fx(2)-fxp(3))^2)+(fx(3)-fxp(5))^2);
%∑ΩŒªΩ«
fz2=atan2((fx(2)-fxp(3)),(fx(1)-fxp(1)));
%∏©—ˆΩ«
fz3=atan2((fx(3)-fxp(5)),sqrt((fx(1)-fxp(1))^2+((fx(2)-fxp(3))^2)));
end