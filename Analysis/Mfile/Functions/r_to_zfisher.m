function  Z = r_to_zfisher(R);
%To convert from correlation to Z-Fisher
Z = 0.5 * log((1+R)./(1-R));
