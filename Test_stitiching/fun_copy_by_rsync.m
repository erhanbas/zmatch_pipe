function exit_code = fun_copy_by_rsync(source_fp, target_fp, verboseQ, copy_folderQ)

if nargin < 3
    verboseQ = false;
    copy_folderQ = false;
elseif nargin < 4
    copy_folderQ = false;
end
if ~copy_folderQ
    if ~endsWith(source_fp, filesep)
        source_fp = [source_fp, filesep];
    end
end
if verboseQ
    cmd_str = sprintf('rsync -rav %s %s', source_fp, target_fp);
    fprintf('Copying %s to %s\n', source_fp, target_fp);
    tmp_tic = tic;
else
    cmd_str = sprintf('rsync -ra %s %s', source_fp, target_fp);
end

system(cmd_str);
if verboseQ
    fprintf('Finish copying. Elapsed time is %s seconds\n', toc(tmp_tic));
end
exit_code = 0;


end