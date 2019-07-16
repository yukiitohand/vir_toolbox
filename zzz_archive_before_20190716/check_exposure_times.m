global LUT_SCTIME_LIST


volumeList = {LUT_SCTIME_LIST.volumeid_IR_1B};


isVesta_I1B = cellfun(@(x) strcmpi(x,'DWNVVIR_I1B_v2'),volumeList);

LUT_SCTIME_LIST_vesta_I1B = LUT_SCTIME_LIST(isVesta_I1B);

expo_times = zeros(length(LUT_SCTIME_LIST_vesta_I1B),1);

for i=1:length(LUT_SCTIME_LIST_vesta_I1B)
    sctimei = LUT_SCTIME_LIST_vesta_I1B(i).sctime;
    vir_obs = get_vir_obs_info(sctimei);
    I1Bdata = VIRdata(vir_obs.basename_ir_1B,vir_obs.dirpath_ir_1B);
    expo_times(i) = I1Bdata.lbl.FRAME_PARAMETER{1}(1);
end

expo_times_vesta = expo_times;

%%
vir_downloader('DWNCLVIR_I1B','dwld',2,'dirskip',0);

isCeresL_I1B = cellfun(@(x) strcmpi(x,'DWNCLVIR_I1B'),volumeList);

LUT_SCTIME_LIST_cl_I1B = LUT_SCTIME_LIST(isCeresL_I1B);

expo_times_cl = zeros(length(LUT_SCTIME_LIST_cl_I1B),1);

for i=1:length(LUT_SCTIME_LIST_cl_I1B)
    sctimei = LUT_SCTIME_LIST_cl_I1B(i).sctime;
    vir_obs = get_vir_obs_info(sctimei);
    I1Bdata = VIRdata(vir_obs.basename_ir_1B,vir_obs.dirpath_ir_1B);
    expo_times_cl(i) = I1Bdata.lbl.FRAME_PARAMETER{1}(1);
end
