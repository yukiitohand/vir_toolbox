function [sctime_info_mp] = get_sctimes_mission_phase(sctime_info,subdir_field,mission_phase_ptr)

is_mp = false(length(sctime_info));
for i=1:length(sctime_info)
    subdir = sctime_info(i).(subdir_field);
    [mission_phase_period] = get_mission_phases(subdir);
    if ~isempty(regexpi(mission_phase_period.phase_name,mission_phase_ptr,'ONCE'))
        is_mp(i) = true;
    end
    
end
sctime_info_mp = sctime_info(is_mp);

end