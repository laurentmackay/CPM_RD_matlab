function out = fortran_subroutine(nm,pars,pars_hdr,local,body)
nm=char(nm);
if nargin==1
    pars='';
    pars_hdr='';
    local={};
    body='';
else
    if ~isempty(pars)
        pars=['(' char(pars) ')'];
    end
    if ~isempty(pars_hdr)
        pars_hdr = [newline char(pars_hdr)];
    end
    if ~isempty(body)
        body=[newline newline char(body) ];
    end
end


if ~isempty(local)
    local_hdr = [newline 'DOUBLE PRECISION ' char(strjoin(local,', ')) ];
else
    local_hdr='';
end

 
out=['SUBROUTINE '  nm  pars  pars_hdr local_hdr body newline 'END SUBROUTINE ' nm newline newline];




end

