function dir = results_dir(expmnt)
global experiment

if nargin==0
    expmnt=experiment;
end

dir = strcat(work_dir(),'results', filesep);

if ~isempty(expmnt)
    dir = strcat(dir,expmnt,filesep);
end

end

