function varargout = pointmatch_task_vessel(filename,from,to,numcores,exitcode)
%% deploys pointmatch function on array task
if nargin==0
    brain = '2017-10-31';tag='';
    runlocal = 0;
    deployment(brain,runlocal);
elseif nargin==1
    disp('Creating *sh file for cluster. File will be here (might take upto 5min to create the file):')
    brain = filename;
    runlocal = 0;
    deployment(brain,runlocal);
    return
elseif nargin==2 % run on local machine
    disp('RUNNING ON LOCAL MACHINE< MIGHT TAKE A WHILE')
    brain = filename;
    runlocal = from;
    deployment(brain,runlocal);
    return
end

if ischar(from);from = eval(from);end
if ischar(to);to = eval(to);end
if ischar(numcores);numcores = eval(numcores);end
if ischar(exitcode);exitcode = eval(exitcode); end
varargout{1} = exitcode;

% first get file list
fid = fopen(filename);
targetlist = textscan(fid,'%s %s %s %s %s %f %f %f %s %f');
fclose(fid);

parpool(numcores)
parfor idx=from:to
    tile1 = targetlist{1}{idx};
    tile2 = targetlist{2}{idx};
    acqusitionfolder1 = targetlist{3}{idx};
    acqusitionfolder2 = targetlist{4}{idx};
    outfold =targetlist{5}{idx};
    pixshift = [targetlist{6}(idx) targetlist{7}(idx) targetlist{8}(idx)];
    ch = targetlist{9}{idx};
    maxnumofdesc = targetlist{10}(idx);
   
    pointmatch(tile1,tile2,acqusitionfolder1,acqusitionfolder2,outfold,pixshift,ch,maxnumofdesc,exitcode)
end
delete(gcp('nocreate'))


end

function deployment(brain,runlocal)
%%
% addpath(genpath('./functions'))
% % addpath(genpath('./common'))
% matlab_cluster_path = '/usr/local/matlab-2017a';
% compiledfunc = '/groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/pointmatch_task/pointmatch_task';
% [mypath,mycomp] = fileparts(compiledfunc);
% taskwrapper = fullfile(pwd,'cluster2.sh');
% taskscript = fullfile(mypath,sprintf('run_%s.sh',mycomp));
clc;clear;
runlocal = true;
brain = '4_17_57_cube';
addpath(genpath('./functions'))
% addpath(genpath('./common'))
matlab_cluster_path = '/usr/local/MATLAB/R2018b';
compiledfunc = '/home/dklab/Documents/Github/Compiled_functions/pointmatch_task/pointmatch_task';
[mypath,mycomp] = fileparts(compiledfunc);
taskwrapper = fullfile(pwd,'cluster2.sh');
taskscript = fullfile(mypath,sprintf('run_%s.sh',mycomp));

% if ~exist(fileparts(compiledfunc),'dir')
%     mkdir(fileparts(compiledfunc));
%     mfilename_ = mfilename('fullpath');
%     unix(sprintf('mcc -m -v -R -I %s -d %s -a %s',mfilename_,fileparts(compiledfunc),fullfile(fileparts(mfilename_),'functions')))
%     %,fullfile(fileparts(mfilename_),'common')
% end

% inputfolder = sprintf('/groups/mousebrainmicro/mousebrainmicro/data/%s/Tiling',brain);
% experimentfolder = sprintf('/nrs/mouselight/cluster/classifierOutputs/%s',brain);
% descriptorfolder = fullfile(experimentfolder,'classifier_output');
% matfolder = fullfile(experimentfolder,'matfiles/');
% scopefile = fullfile(matfolder,'scopeloc.mat');

inputfolder = sprintf('/data/Vessel/ML_stitching/%s/raw_data',brain);
experimentfolder = sprintf('/data/Vessel/ML_stitching/%s',brain);
descriptorfolder = fullfile(experimentfolder,'stage_2_descriptor_output');
matfolder = fullfile(experimentfolder,'matfiles/');
scopefile = fullfile(matfolder,'scopeloc.mat');

if ~exist(matfolder,'dir')
    mkdir(matfolder)
end
if exist(scopefile,'file')
    load(scopefile,'scopeloc')
    [neighbors] = buildNeighbor(scopeloc.gridix(:,1:3));
else
    newdash = 1; % set this to 1 for datasets acquired after 160404
    [scopeloc] = getScopeCoordinates(inputfolder,newdash);% parse from acqusition files
    [neighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
    imsize_um = [384.72339, 456.35, 250];
    save(scopefile,'neighbors','scopeloc','imsize_um','experimentfolder','inputfolder')
end
directionMap = containers.Map({'-X','-Y','X','Y','-Z','Z'},[ 2, 3, 4, 5, 6, 7]);
%% Load feature match in Z direction and test if unreliable tiles exist. 
directions = 'Z';
ch='0';
maxnumofdesc = 10e3;
checkversion=1;
if 0
    pixinit = zeros(size(neighbors,1),3);
    nummatches = zeros(size(neighbors,1),1);
else
    % load finished tile matches. find badly matched or missing tile pairs
    [regpts,featmap] = loadMatchedFeatures(scopeloc,descriptorfolder,directions,checkversion);
    % initalize missing tiles based on knn
    numthr = 50;
    [pixinit,nummatches] = initTiles(featmap,directions,scopeloc,numthr);
    fprintf('Total number of matches in Z direction: %d\n', sum(nummatches));
end
% badtiles is a logical column vector. True if tile has bad matching in the
% z positive direction
badtiles = nummatches<numthr & ~isnan(neighbors(:,directionMap(directions)));
%%
% generate a txt file with all tile pairs to be matched with optional shift
% values
% (filename,from,to,numcores,exitcode)
% pixinit is empty from the steps above. 
if ~exist('pixinit', 'var') || isempty(pixinit)
    pixinit = zeros(size(neighbors,1),3);
end
outlistfile = fullfile(pwd,'shfiles',sprintf('outlistfile_%s_%s.txt',brain,date));
if ~runlocal;fid = fopen(outlistfile,'w');end
for ii = 1:length(badtiles)
    if ~badtiles(ii)
        continue
    end
    tile1 = fullfile(descriptorfolder,scopeloc.relativepaths{ii});
    acqusitionfolder1 = fileparts(scopeloc.filepath{ii});
    iineig = neighbors(ii,directionMap(directions));
    tile2 = fullfile(descriptorfolder,scopeloc.relativepaths{iineig});
    acqusitionfolder2 = fileparts(scopeloc.filepath{iineig});
    outfold = tile1;
    if runlocal
        pointmatch_vessel(tile1, tile2, acqusitionfolder1, acqusitionfolder2, outfold, pixinit(ii,:), ch, maxnumofdesc,0);
    else
        fprintf(fid,'%s %s %s %s %s %f %f %f %s %f\n',tile1,tile2,acqusitionfolder1,acqusitionfolder2,outfold,pixinit(ii,:),ch,maxnumofdesc);
    end
end
if ~runlocal;fclose(fid);end
%%

if ~runlocal
    %% task script (filename,from,to,numcores,exitcode)
    %find number of random characters to choose from and %specify length of random string to generate
    s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';numRands = length(s);sLength = 10;
    %-o /dev/null
    esttime = 30*60;
    numcores = 16;
    exitcode = 0;
    tottasks = sum(badtiles(:));
    intervals = round(linspace(0,tottasks,round(tottasks/numcores)));
    
    mkdir(fullfile(pwd,'shfiles'))
    myfile = fullfile(pwd,'shfiles',sprintf('featmatchrun_%s_%s.sh',brain,date))
    if ~runlocal; fid = fopen(myfile,'w'); end
    for ii = 1:length(intervals)-1
        from = intervals(ii)+1;
        to = intervals(ii+1);
        cmd = sprintf('''%s %s %s %s %d %d %d %d''',taskwrapper,taskscript,matlab_cluster_path,outlistfile,from,to,numcores,exitcode);
        
        %generate random string
        randString = s( ceil(rand(1,sLength)*numRands) );
        name = sprintf('zm_%05d-%s',ii,randString);
        mysub = sprintf('bsub -J %s -n%d -R"affinity[core(1)]" -We %d -o /dev/null %s\n',name,numcores,esttime/60,cmd);
        fwrite(fid,mysub);
    end
    if ~runlocal; fclose(fid); unix(sprintf('chmod +x %s',myfile)); disp(myfile);end
end
end




