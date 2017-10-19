function featmap = loadMatchedFeatures(descriptorfolder,directions,featmap)
% get list file
myfile = dir(fullfile(descriptorfolder,'list*files'));
fid=fopen(fullfile(descriptorfolder,myfile.name),'r');
inputfiles = textscan(fid,'%s');inputfiles = inputfiles{1};fclose(fid);
folders = unique(cellfun(@fileparts,inputfiles,'UniformOutput',false),'stable');

if nargin<2
    directions='Z'; % only check z
    for idx=1:length(folders)
        [featmap(idx).X,featmap(idx).Y,featmap(idx).Z] = deal([]);
    end
elseif nargin<3
    for idx=1:length(folders)
        [featmap(idx).X,featmap(idx).Y,featmap(idx).Z] = deal([]);
    end
end
%%
tmp=cell(1,length(folders));
parfor idx = 1:length(folders)
    matchfile = fullfile(folders{idx},sprintf('match-%s.mat',directions));
    if exist(matchfile,'file')
        % load descriptors
        tmp{idx} = load(matchfile)
    end
end
%%
for idx = 1:length(folders)
    featmap(idx).(genvarname(directions)) = tmp{idx};
end
