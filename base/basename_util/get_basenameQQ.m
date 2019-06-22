function [basenameQQ] = get_basenameQQ(basenameQUB)
    prop = getProp_basenameVIR(basenameQUB);
    if ~isempty(prop)
        prop_HK = prop;
        prop_HK.type = 'QQ';
        basenameQQ = get_basenameVIR_fromProp(prop_HK);
    else
        basenameQQ = '';
    end

end