function [basenameHK] = get_basenameHK(basenameQUB)
    prop = getProp_basenameVIR(basenameQUB);
    if ~isempty(prop)
        prop_HK = prop;
        prop_HK.type = 'HK';
        basenameHK = get_basenameVIR_fromProp(prop_HK);
    else
        basenameHK = '';
    end

end