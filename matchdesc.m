function matchdesc(descriptorfile,scopefile,outfolder,indstart,indend,thr,numcores,zshift)
%ZMATCH Summary of this function goes here
%
% [OUTPUTARGS] = ZMATCH(INPUTARGS) Explain usage here
%
% Inputs:
%
% Outputs:
%
% Examples:
%
% Provide sample usage code here
%
% See also: List related files here

% $Author: base $	$Date: 2016/09/30 18:51:23 $	$Revision: 0.1 $
% Copyright: HHMI 2016
if nargin<1
    brain = '2017-09-19';
    tag='';
    runlocal = 0;
    deployment(brain,tag,runlocal)
    return
end
if ~isdeployed
    addpath(genpath('./thirdparty'))
    addpath(genpath('./functions'))
end
indstart = str2double(indstart);
indend = str2double(indend);
thr = str2double(thr);
numcores = str2double(numcores);
zshift = str2double(zshift);
% addpath(genpath('./thirdparty'))
outfile = fullfile(outfolder,sprintf('%05d_%05d-pointmatch',indstart,indend));
% load descriptor file
load(descriptorfile,'descriptors')
% load scopelocation file
load(scopefile,'scopeloc','imsize_um','experimentfolder','inputfolder')
[neighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
checkthese = [1 4 5 7]; % 0 - below
neigs = neighbors(indstart:indend,checkthese);

paireddescriptor = arraytask(descriptors,neigs,scopeloc,thr,numcores,zshift);

save(outfile,'paireddescriptor')
end
function deployment(brain,tag,runlocal)
% 
compiledfunc = '/groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/matchdesc/matchdesc'
outfold = fileparts(compiledfunc);
mkdir(outfold)
if 0
%     mcc -m -v -R -singleCompThread /groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/zmatch_pipe/matchdesc.m -d /groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/matchdesc -a /groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/zmatch_pipe/functions
%     mcc -m -v -R -I /groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/zmatch_pipe/matchdesc.m -d /groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/matchdesc_mult -a /groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/zmatch_pipe/functions
    unix(sprintf('mcc -m -v -R -singleCompThread /groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/zmatch_pipe/matchdesc.m -d %s -a /groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/zmatch_pipe/functions',outfold))
    sprintf('mcc -m -v -R -I /groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/zmatch_pipe/matchdesc.m -d %s -a /groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/zmatch_pipe/functions',outfold)
end
%%
% brain = '2015-06-19';
% tag = '_backup'
if nargin<3
    runlocal = 1;
end
experimentfolder = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/%s%s/',brain,tag)
matfolder = fullfile(experimentfolder,'matfiles/');
descriptorfile = fullfile(matfolder,'descriptors_ch0');
scopefile = fullfile(matfolder,'scopeloc');
% scopeparams = fullfile(matfolder,'scopeparams');
% scopeparams = fullfile(matfolder,'scopeparams_pertile');

outfolder = fullfile(matfolder,'pointmatches2')
mkdir(outfolder)

numcores = 8;
mkdir(fullfile(pwd,'shfiles'))
myfile = fullfile(pwd,'shfiles',sprintf('zmatchrun_%s_%s.sh',brain,date))

%find number of random characters to choose from
s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
numRands = length(s);
%specify length of random string to generate
sLength = 10;
%-o /dev/null
esttime = 30*60;

load(fullfile(matfolder,'scopeloc'),'scopeloc')
Ntiles = size(scopeloc.loc,1);

% inds = round(linspace(0,Ntiles,Ntiles/100+1));
inds = round(linspace(0,Ntiles,1000+1));
inds = [0:numcores:Ntiles-1 Ntiles];
% inds = 729:855;%round(linspace(729,855,127));

if 1
    myfiles = dir(fullfile(outfolder,'*.mat'));
    doneinds = cellfun(@(x) str2num(x(1:5)),{myfiles.name});
    [finished,bb] = min(pdist2((inds+1)',doneinds(:)),[],2);finished = ~finished;
    if isempty(finished)
        finished = zeros(1,length(inds)-1);
    end
    % [finished,bb] = min(pdist2((inds+1)',find(cellfun(@isempty,regpts))'),[],2)
else
    finished = zeros(1,length(inds)-1);
    %     finished(402:494) = 1;
    %     finished = ~finished
end
sum(~finished)

%%

% runlocal=0
if ~runlocal
    fid = fopen(myfile,'w');
end
% zmatch(descriptorfile,scopefile,scopeparams,outfolder,indstart,indend)
% for ii=1:length(inds)
thr = .1
zshift = 0
for ii=1:length(inds)-1
    %%
    if finished(ii)
        continue
    end
    ii
    %generate random string
    randString = s( ceil(rand(1,sLength)*numRands) );
    name = sprintf('zm_%05d-%s',ii,randString);
    args = sprintf('''%s %s %s %s %05d %05d %0.1f %d %d''',compiledfunc,descriptorfile,scopefile,outfolder,inds(ii)+1,inds(ii+1),thr,numcores,zshift);
%     mysub = sprintf('qsub -pe batch %d -l d_rt=%d -N %s -j y -o /dev/null -b y -cwd -V %s\n',numcores,esttime,name,args);
    mysub = sprintf('bsub -n%d -We %d -J %s -o /dev/null %s\n',numcores,esttime/60,name,args);
    if runlocal
        matchdesc(sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/%s/matfiles/descriptors_ch0',[brain,tag]),...
            sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/%s/matfiles/scopeloc',[brain,tag]),...
            sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/%s/matfiles/pointmatches',[brain,tag]),...
            sprintf('%05d',inds(ii)+1), sprintf('%05d',inds(ii+1)), '0.1',num2str(feature('numCores')),'0')
    else
        fwrite(fid,mysub);
    end
end
unix(sprintf('chmod +x %s',myfile));
if ~runlocal
    fclose(fid);
    myfile
end
myfile

end
function paireddescriptor=arraytask(descriptors,neigs,scopeloc,thr,numcores,zshift)
dims = [1024 1536 251]; % in xyz
paireddescriptor = cell(size(neigs,1),1);
maxnumofdesc=3e3;
projectionThr = 5;

%%
tag = 'XYZ';
% scopefile2 = readScopeFile(acqusitionfolder2);
% imsize_um = [scopefile1.x_size_um,scopefile1.y_size_um,scopefile1.z_size_um];
% % estimate translation
% gridshift = ([scopefile2.x scopefile2.y scopefile2.z]-[scopefile1.x scopefile1.y scopefile1.z]);
% iadj =find(gridshift);
% stgshift = 1000*([scopefile2.x_mm scopefile2.y_mm scopefile2.z_mm]-[scopefile1.x_mm scopefile1.y_mm scopefile1.z_mm]);
% pixshift = round(stgshift.*(dims-1)./imsize_um);

%%
delete(gcp('nocreate'))
parpool(numcores)
%%
parfor ineig = 1:size(neigs,1)
    paireddescriptor{ineig}.X = [];
    paireddescriptor{ineig}.Y = [];
    paireddescriptor{ineig}.neigs = neigs(ineig,:);
    paireddescriptor{ineig}.matchrate = 0;
end
%%
parfor ineig = 1:size(neigs,1)
    %% load descriptor pairs X (center) - Y (adjacent tile)
    idxcent = neigs(ineig,1);
    idxadj =  neigs(ineig,end);
    if isnan(idxadj)
        continue
    end
    
    scopefile1 = readScopeFile(fileparts(scopeloc.filepath{idxcent}));
    scopefile2 = readScopeFile(fileparts(scopeloc.filepath{idxadj}));
    imsize_um = [scopefile1.x_size_um,scopefile1.y_size_um,scopefile1.z_size_um];
    % estimate translation
    gridshift = ([scopefile2.x scopefile2.y scopefile2.z]-[scopefile1.x scopefile1.y scopefile1.z]);
    iadj =find(gridshift);
    stgshift = 1000*([scopefile2.x_mm scopefile2.y_mm scopefile2.z_mm]-[scopefile1.x_mm scopefile1.y_mm scopefile1.z_mm]);
    pixshift = round(stgshift.*(dims-1)./imsize_um);
    
    % read descs
    %     desc1 = readDesc(tile1);
    %     desc2 = readDesc(tile2);
    desc1 = descriptors{idxcent};
    desc2 = descriptors{idxadj};
    % check if input exists
    if isempty(desc1) | isempty(desc2)
        rate_ = 0;
        X_ = [];
        Y_ = [];
        uni = 0;
    else
        % correct images, xy flip
        desc1 = correctTiles(desc1,dims);
        desc2 = correctTiles(desc2,dims);
        % truncate descriptors
        desc1 = truncateDesc(desc1,maxnumofdesc,thr);
        desc2 = truncateDesc(desc2,maxnumofdesc,thr);
        % idaj : 1=right(+x), 2=bottom(+y), 3=below(+z)
        % pixshift(iadj) = pixshift(iadj)+expensionshift(iadj); % initialize with a relative shift to improve CDP
        clc
        debug = 0;
        matchparams = modelParams(projectionThr,debug);
        %%
        if length(iadj)~=1 | max(iadj)>3
            error('not 6 direction neighbor')
        end
        [X_,Y_,out,rate_,pixshiftout,nonuniformity] = searchpair(desc1(:,1:3),desc2(:,1:3),pixshift,iadj,dims,matchparams);
        if ~isempty(X_)
            X_ = correctTiles(X_,dims);
            Y_ = correctTiles(Y_,dims);
        else % try a more conservative run
            matchparams_ = matchparams;
            matchparams_.opt.outliers = .5;
            [X_,Y_,out,rate_,pixshiftout,nonuniformity] = searchpair(desc1(:,1:3),desc2(:,1:3),pixshift,iadj,dims,matchparams_);
            if ~isempty(X_)
                X_ = correctTiles(X_,dims);
                Y_ = correctTiles(Y_,dims);
            end
        end
        uni = mean(nonuniformity)<=.5;
    end
    paireddescriptor{ineig}.matchrate = rate_;
    paireddescriptor{ineig}.X = X_;
    paireddescriptor{ineig}.Y = Y_;
    paireddescriptor{ineig}.uni = uni;
end
delete(gcp('nocreate'))
end
function [neighbors] = buildNeighbor(grididx)
%BUILDNEIGHBOR Given grid location, builds an connectivity matrix and
%neighbor edge graph
%
% [NEIGHBORS] = BUILDNEIGHBOR(GRIDIDX)
%
% Inputs: [grididx]: Nx3: tile subscripts locations
%
% Outputs: [neighbors]: Nx7: neighbors list in [idx left top right bottom
% above below]:[id -x -y +x +y -z +z] format
%
% Examples:
% [aa,bb,cc] = ndgrid([-1:1],[-1:1],[-1:1]);
% grididx = 5+[aa(:),bb(:),cc(:)];
% grididx(1:3:end,:) = [];
% [neighbors] = buildNeighbor(grididx)
%
% See also:

% $Author: base $	$Date: 2016/08/19 11:01:15 $	$Revision: 0.1 $
% Copyright: HHMI 2016

dims = max(grididx);
N = size(grididx,1);
neighbors = nan(N,7); %[idx left top right bottom above below]:[id -x -y +x +y -z +z]
sixneig = [
    [-1 0 0]; % left (-x)
    [0 -1 0]; % top (-y)
    [1 0 0]; % right (+x)
    [0 1 0]; % bottom (+y)
    [0 0 -1]; % above (-z)
    [0 0 1]]; % below (+z)
indsgrid = sub2ind(dims,grididx(:,1),grididx(:,2),grididx(:,3));
indIM = NaN(dims);
indIM(indsgrid) = 1:N;
for idxN = 1:N
    sub = grididx(idxN,:);
    % check 6 neighbor
    set = nan(6,1);
    sub_ = nan(6,3);
    for ix = 1:6
        sub_(ix,:) = sub + sixneig(ix,:);
    end
    validsubs = all(sub_>0&sub_<=ones(6,1)*dims,2);
    sub__ = sub_(validsubs,:);
    set(validsubs) = indIM(sub2ind(dims,sub__(:,1),sub__(:,2),sub__(:,3)));
    neighbors(idxN,:) = [idxN;set];
end
end









