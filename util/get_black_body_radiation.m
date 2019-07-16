function [bb_spc] = get_black_body_radiation(wv_um,T)
% [bb_spc] = get_black_body_radiation(wv,T)
%  INPUTS
%   wv_um: wavelength vector, scalar [micron]
%   T : temperature Kelvin
%  OUTPUTS
%   bb_spc: black body radiation spectrum [W / (sr*m^2*m )]

% planck constant
h = 6.62607015*10^(-34); % Js

% speed of light
c = 299792458; % m/s

% Boltzman constant
kB = 1.380649*10^(-23); % J/K

% frequency of light
wv_m = wv_um*10^(-6); %
% f = c * (10^6) ./ wv; % 1/s

bb_spc = (2*h*c^2)./ (wv_m.^5) ./ (exp((h*c)./(kB*T.*wv_m)) - 1 );

end