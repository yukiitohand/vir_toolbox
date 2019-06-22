function [lnks] = get_links_remoteHTML_vir(html)
% [lnks] = get_links_remoteHTML_vir(html)
%   match files in the html files obtained at the remote server 
% Input
%   html: html text
% Output
%   lnks: struct having two fields
%     hyperlink: hyperlink, only the basename (or with extention) is read.
%                no slash will be attached
%     type     : type of the content at the link {'PARENTDIR','To Parent Directory','dir'}
%                if other types are specified, it will be regarded as a
%                file.

ptrn_lnk = '<img src="[^<>"]*" alt="\[(?<type>[^<>"]*)\]"></td><td><a href="(?<hyperlink>[^<>"]+)">';
lnks = regexpi(html,ptrn_lnk,'names');

for i=1:length(lnks)
    if strcmpi(lnks(i).type,'PARENTDIR')
        lnks(i);
    else
        lnk = strip(lnks(i).hyperlink,'right','/');
        [~,link_name,ext] = fileparts(lnk);
        link_name2 = [link_name ext];
        lnks(i).hyperlink = link_name2;
    end
end

end