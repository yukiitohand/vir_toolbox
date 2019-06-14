function [hdr_info] = extract_imghdr_from_lbl_vir(lbl_info)
% [hdr_info] = extract_imghdr_from_lbl_vir(lbl_info)
%   extract header information (envi format) from VIR LABEL file
%  Input Parameters
%   lbl: struct of LABEL file
%  Output Parameters
%   hdr_info: struct of header in envi format, if no image is found, [] is
%             returend.

obj_file_image = lbl_info.OBJECT_QUBE;

if isempty(obj_file_image)
    hdr_info = [];
else
    hdr_info = [];
    hdr_info.samples = obj_file_image.CORE_ITEMS{1}(2);
    hdr_info.lines = obj_file_image.CORE_ITEMS{1}(3);
    hdr_info.bands = obj_file_image.CORE_ITEMS{1}(1);

    if strcmpi(obj_file_image.CORE_ITEM_TYPE,'IEEE_REAL')
        hdr_info.data_type = 4;
        hdr_info.byte_order = 1;
    elseif strcmp(obj_file_image.CORE_ITEM_TYPE,'MSB_INTEGER')
        hdr_info.byte_order = 1;
        if obj_file_image.CORE_ITEM_BYTES==2
            hdr_info.data_type = 2;
        elseif obj_file_image.CORE_ITEM_BYTESS==1
            % maybe?
            hdr_info.data_type = 1;
        else
            error('Undefined "obj_file_image.CORE_ITEM_BYTES"');
        end
    else
        error('The data type: %s is not supported.',obj_file_image.CORE_ITEM_TYPE);
    end

    hdr_info.header_offset = 0;
    % hdr_info.header_offset = img_obj.RECORD_BYTES;
    hdr_info.wavelength = obj_file_image.BAND_BIN_CENTER{1}(:);
    hdr_info.fwhm = obj_file_image.BAND_BIN_WIDTH{1}(:);

    hdr_info.interleave = 'bip';

end




end