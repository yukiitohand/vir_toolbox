function [I1Bim_if] = vir_ra2if(I1Bdata,SolSpcdata)


if isempty(I1Bdata.img)
    I1Bdata.readimg();
end

d_km = I1Bdata.lbl.SPACECRAFT_SOLAR_DISTANCE{1};
[ d_au ] = km2au( d_km );
if isempty(SolSpcdata.tab)
    SolSpcdata.readTAB();
end

SS= [SolSpcdata.tab.data.SPECTRAL_IRRADIANCE]';
B = length(SS);

I1Bim_if = I1Bdata.img  .* ((pi .* (d_au.^2)) ./ reshape(SS,[1,1,B]));

end