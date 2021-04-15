function set_protocol(nm)
global RD_base protocol

if isempty(RD_base)
    get_RD_base();
end

if exist(strcat(RD_base,'protocols',filesep, nm),'file')
    if ~isempty(protocol)
        rmpath(strcat(RD_base,'protocols',filesep, nm));
    end
    protocol=nm;
    addpath(strcat(RD_base,'protocols',filesep, nm));
else
    error(strcat('Could not find "',nm,'" in the protocols folder.'));
end

end

