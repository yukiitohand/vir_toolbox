function [mission_phase_period] = get_mission_phases(subdir)
% [mission_phase_period] = get_mission_phases(subdir)
%  get sub information for mission from subdirectory something like
%      'DATA/20110503_APPROACH/20110510_OPNAV_002'
%  INPUT: subdir: path
%  OUTPUT: mission_phase_period, struct having fields:
%    phase_date
%    phase_name
%    period_date
%    period_name

phase_ptrn = 'DATA/(?<phase_date>[\d]{8})_(?<phase_name>[^/]*)/*(?<period_date>[\d]{8})*_*(?<period_name>.*)*';

mission_phase_period = regexpi(subdir,phase_ptrn,'names');

end