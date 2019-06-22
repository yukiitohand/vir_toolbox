classdef VIRdata < HSI
    %CRISMdata class
    %   For any type of CRISM data (EDR,CDR,TRDR,TER,...)
    
    properties
        lblpath;
        tabpath;
        lbl;
        tab;
        hkt;
        qq;
        prop = [];
        data_type = '';
    end
    
    methods
        function obj = VIRdata(basename,dirpath,varargin)
            % load property and find out the type of data from "basename"
            prop = getProp_basenameVIR(basename);
            if ~isempty(prop)
                if isempty(prop.type)
                    data_type = 'QUBE';
                elseif strcmpi(prop.type,'HK')
                    data_type = 'HK';
                elseif strcmpi(prop.type,'QQ')
                    data_type = 'QQ';
                else
                    error('Undefined data type');
                end
            else
                prop = getProp_basenameVIRCAL(basename);
                data_type = prop.name;
            end
                
            % find out extension

            obj@HSI(basename,dirpath);
            obj.lblpath = joinPath(dirpath,[basename '.LBL']);

            obj.lbl = crismlblread_v2(obj.lblpath);
            
            switch data_type
                case {'QUBE'}
                    obj.hdr = extract_imghdr_from_lbl_vir(obj.lbl);
                    obj.imgpath = joinPath(dirpath,[basename '.QUB']);
                case {'QQ'}
                    obj.hdr = extract_qqhdr_from_lbl_vir(obj.lbl);
                    obj.imgpath = joinPath(dirpath,[basename '.QUB']);
                case {'HK'}
                    obj.hdr = [];
                    obj.tabpath = joinPath(dirpath,[basename '.TAB']);
                case {'RESP'}
                    obj.hdr = extract_dathdr_from_lbl_vir(obj.lbl);
                    obj.imgpath = joinPath(dirpath,[basename '.DAT']);
                otherwise
                    obj.hdr = [];
                    obj.tabpath = joinPath(dirpath,[basename '.TAB']);
            end
            
            obj.data_type = data_type;
            obj.prop = prop;

        end
        function [tab] = readTAB(obj)
            [ tab ] = crismTABread( obj.dirpath, obj.lbl );
            if isempty(tab)
                fprintf('no tab is found');
            end
            obj.tab = tab;
        end
        function [] = loadHKTAB(obj)
            basenameHK = get_basenameHK(obj.basename);
            obj.hkt = VIRdata(basenameHK,obj.dirpath);
        end
        function [] = loadQQ(obj)
            basenameQQ = get_basenameQQ(obj.basename);
            obj.qq = VIRdata(basenameQQ,obj.dirpath);
        end
    end
end