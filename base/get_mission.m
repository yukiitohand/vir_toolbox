function [mission] = get_mission(mission_id)
% [mission] = get_mission(mission_id)
%  get mission target from mission ID
%   INPUTS: mission_id {'C','V','X'}
%   OUTPUTS: 
%     'Ceres'  for 'C'
%     'Vesta'  for 'V'
%     'Cruise' for 'X'

switch upper(mission_id)
    case 'C'
        mission = 'Ceres';
    case 'V'
        mission = 'Vesta';
    case 'X'
        mission = 'Cruise';
    otherwise
        error('Mission is right??');
end

end