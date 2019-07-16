function [dark_s] = stretch_dark(dark,s)

dark_mean = nanmean(dark);
dark_meansub = dark - dark_mean;
dark_s = dark_mean + s*dark_meansub;

end