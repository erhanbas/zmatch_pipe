function varargout = pointmatch_array(filename,from,to,numcores,exitcode)
%% deploys pointmatch function on array task
if nargin<1
    brain = '2017-09-25';tag='';runlocal = 0;deployment(brain,tag,runlocal);return
end
% if ischar(pixshift)
%     pixshift = eval(pixshift); % pass initialization
% end

if 1
    % first get file list
    fid = fopen(filename); 
    targetlist = textscan(fid,'%s %s %d %d %d','Delimiter','\n');
    targetlist=targetlist{1}; 
    fclose(fid);

else
    % pointmatch_array(descriptorfolder,scopefile,from,to,pixshift,numcores,exitcode)
    tr=load(scopefile,'scopeloc');scopeloc=tr.scopeloc;
    [neighbors] = buildNeighbor(scopeloc.gridix(:,1:3));
    delete(gcp('nocreate'))
    parpool(numcores)
    parfor idx=from:to
        zidx = neighbors(idx,end);
        if isnan(zidx);continue;end
        tile1 = fullfile(descriptorfolder,scopeloc.relativepaths{idx});
        tile2 = fullfile(descriptorfolder,scopeloc.relativepaths{zidx});
        acqusitionfolder1 = fileparts(scopeloc.filepath{idx});
        acqusitionfolder2 = fileparts(scopeloc.filepath{zidx});
        outfold =tile1;
        % initialize pix based on preruns
        % initpix()
        % pixshift ='[0 0 0]';
        pointmatch(tile1,tile2,acqusitionfolder1,acqusitionfolder2,outfold,pixshift,exitcode)
    end
    delete(gcp('nocreate'))
end


end

function deployment(brain,tag,runlocal)
%%
addpath(genpath('./functions'))
addpath(genpath('./common'))

inputfolder = sprintf('/groups/mousebrainmicro/mousebrainmicro/data/%s/Tiling',brain);
experimentfolder = sprintf('/nrs/mouselight/cluster/classifierOutputs/%s',brain);

descriptorfolder = fullfile(experimentfolder,'classifier_output');
matfolder = fullfile(experimentfolder,'matfiles/');
scopefile = fullfile(matfolder,'scopeloc.mat');

if ~exist(matfolder,'dir')
    mkdir(matfolder)
end
if exist(scopefile,'file')
    load(scopefile,'scopeloc')
else
    newdash = 1; % set this to 1 for datasets acquired after 160404
    [scopeloc] = getScopeCoordinates(inputfolder,newdash);% parse from acqusition files
    save(fullfile(matfolder,'scopeloc'),'scopeloc','imsize_um','experimentfolder','inputfolder')
end
[neighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
directionMap = containers.Map({'-X','-Y','X','Y','-Z','Z'},[ 2, 3, 4, 5, 6, 7]);
%%
directions = 'Z';
if 0
    pixinit = zeros(size(neighbors,1),3);
    nummatches = zeros(size(neighbors,1),1);
else
    % load finished tile matches. find badly matched or missing tile pairs
    featmap = loadMatchedFeatures(descriptorfolder,directions);
    % initalize missing tiles based on knn
    numthr = 50;
    [pixinit,nummatches] = initTiles(featmap,directions,scopeloc,numthr);
end

%%
% generate a txt file with all tile pairs to be matched with optional shift
% values
badtiles = nummatches<numthr;
% (filename,from,to,numcores,exitcode)
outlistfile = 'mytest.txt';
fid = fopen(outlistfile,'w');
for ii = find(badtiles(:)')
    tile1 = fullfile(descriptorfolder,scopeloc.relativepaths{ii});
    iineig = neighbors(ii,directionMap(directions));
    if ~isnan(iineig)
        tile2 = fullfile(descriptorfolder,scopeloc.relativepaths{iineig});
        fprintf(fid,'%s %s %d %d %d\n',tile1,tile2,pixinit(ii,:));
    end
end
fclose(fid);

%% task script

%%
% /groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/pointmatch_array

numcores = 16;
exitcode = 0;
mkdir(fullfile(pwd,'shfiles'))
myfile = fullfile(pwd,'shfiles',sprintf('featmatchrun_%s_%s.sh',brain,date))
compiledfunc = '/groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/pointmatch_array/pointmatch_array'
if ~exist(fileparts(compiledfunc),'dir');
    mkdir(fileparts(compiledfunc));
    mfilename_ = mfilename('fullpath');
    unix(sprintf('mcc -m -v -R -I %s -d %s -a %s',mfilename_,fileparts(compiledfunc),fullfile(fileparts(mfilename_),'functions')))
end

%find number of random characters to choose from and %specify length of random string to generate
s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';numRands = length(s);sLength = 10;
%-o /dev/null
esttime = 30*60;


%%
Ntiles = size(scopeloc.loc,1);
% inds = round(linspace(0,Ntiles,1000+1));
inds = [0:numcores:Ntiles-1 Ntiles];
finished = checkmissing(outfolder,inds);
if isempty(finished); finished = zeros(1,length(inds)-1); end
sum(~finished)

%%
% runlocal=0
if ~runlocal; fid = fopen(myfile,'w'); end

for ii=1:length(inds)-1
    %%
    if finished(ii)
        continue
    end
    disp(ii)

    %generate random string
    randString = s( ceil(rand(1,sLength)*numRands) );
    name = sprintf('zm_%05d-%s',ii,randString);
    
    from = inds(ii)+1;
    to = inds(ii+1);
    pixshift = [0 0 0];
    
    args = sprintf('''%s %s %s %d %d "[%d %d %d]" %d %d''',compiledfunc,descriptorfolder,scopefile,from,to,pixshift,numcores,exitcode);
    mysub = sprintf('bsub -n%d -R"affinity[core(1)]" -We %d -J %s -o /dev/null %s\n',numcores,esttime/60,name,args);
    if runlocal
        pointmatch_array(descriptorfolder,scopefile,from,to,pixshift,numcores,exitcode)
    else
        fwrite(fid,mysub);
    end
end

if ~runlocal; fclose(fid); unix(sprintf('chmod +x %s',myfile)); disp(myfile);end

end

function finished = checkmissing(outputfolder,inds)
myfiles = dir(fullfile(outputfolder,'*.mat'));
doneinds = cellfun(@(x) str2num(x(1:5)),{myfiles.name});
[finished,bb] = min(pdist2((inds+1)',doneinds(:)),[],2);finished = ~finished;
end


function deployment_(brain,tag,runlocal)
%%
experimentfolder = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/%s%s/',brain,tag)
matfolder = fullfile(experimentfolder,'matfiles/');
descriptorfolder = sprintf('/nrs/mouselight/cluster/classifierOutputs/%s/classifier_output',brain)
scopefile = fullfile(matfolder,'scopeloc');
load(fullfile(matfolder,'scopeloc'),'scopeloc')
outfolder = fullfile(matfolder,'pointmatches');
mkdir(outfolder)

numcores = 16;
exitcode = 0;
mkdir(fullfile(pwd,'shfiles'))
myfile = fullfile(pwd,'shfiles',sprintf('featmatchrun_%s_%s.sh',brain,date))
compiledfunc = '/groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/pointmatch_array/pointmatch_array'
if ~exist(fileparts(compiledfunc),'dir');
    mkdir(fileparts(compiledfunc));
    mfilename_ = mfilename('fullpath');
    unix(sprintf('mcc -m -v -R -I %s -d %s -a %s',mfilename_,fileparts(compiledfunc),fullfile(fileparts(mfilename_),'functions')))
end

%find number of random characters to choose from and %specify length of random string to generate
s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';numRands = length(s);sLength = 10;
%-o /dev/null
esttime = 30*60;

Ntiles = size(scopeloc.loc,1);
% inds = round(linspace(0,Ntiles,1000+1));
inds = [0:numcores:Ntiles-1 Ntiles];

finished = checkmissing(outfolder,inds);
if isempty(finished); finished = zeros(1,length(inds)-1); end
sum(~finished)

%%
% runlocal=0
if ~runlocal; fid = fopen(myfile,'w'); end

for ii=1:length(inds)-1
    %%
    if finished(ii)
        continue
    end
    disp(ii)

    %generate random string
    randString = s( ceil(rand(1,sLength)*numRands) );
    name = sprintf('zm_%05d-%s',ii,randString);
    
    from = inds(ii)+1;
    to = inds(ii+1);
    pixshift = [0 0 0];
    
    args = sprintf('''%s %s %s %d %d "[%d %d %d]" %d %d''',compiledfunc,descriptorfolder,scopefile,from,to,pixshift,numcores,exitcode);
    mysub = sprintf('bsub -n%d -R"affinity[core(1)]" -We %d -J %s -o /dev/null %s\n',numcores,esttime/60,name,args);
    if runlocal
        pointmatch_array(descriptorfolder,scopefile,from,to,pixshift,numcores,exitcode)
    else
        fwrite(fid,mysub);
    end
end

if ~runlocal; fclose(fid); unix(sprintf('chmod +x %s',myfile)); disp(myfile);end

end




