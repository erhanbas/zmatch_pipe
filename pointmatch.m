function varargout = pointmatch(tile1,tile2,acqusitionfolder1,acqusitionfolder2,outfold,exitcode)
%%
if ~isdeployed
    addpath(genpath('./thirdparty'))
    addpath(genpath('./functions'))
    rmpath(genpath('./common'))
end

if nargin<1
    rawfolder = '/groups/mousebrainmicro/mousebrainmicro/data/'
    classifierfolder = '/nrs/mouselight/cluster/classifierOutputs/'
    sample = '2017-08-28'
    tileid1 = '/2017-08-31/00/00740'
    tileid2 = '/2017-08-31/00/00986'
    tile1 = fullfile(classifierfolder,sample,'/classifier_output',tileid1);
    tile2 = fullfile(classifierfolder,sample,'/classifier_output',tileid2);
    acqusitionfolder1 = fullfile(rawfolder,sample,'Tiling',tileid1);
    acqusitionfolder2 = fullfile(rawfolder,sample,'Tiling',tileid2);
    [neighbors] = buildNeighbor(scopeloc.gridix(:,1:3));
    outfile = fullfile(classifierfolder,sample,'/classifier_output',tileid1)
    imsize_um = [386.67 423.72 250]
end

if nargin<5
    outfold = tile1;
    exitcode = 0;
elseif nargin < 6
    exitcode = 0;
end
varargout{1} = exitcode;
dims = [1024,1536,251];
%%
% read descs
desc1 = readDesc(tile1);
desc2 = readDesc(tile2);
% correct images, xy flip
desc1 = correctTiles(desc1,dims);
desc2 = correctTiles(desc2,dims);
scopefile1 = readScopeFile(acqusitionfolder1);
scopefile2 = readScopeFile(acqusitionfolder2);
imsize_um = [scopefile1.x_size_um,scopefile1.y_size_um,scopefile1.z_size_um];
% estimate translation
gridshift = ([scopefile2.x scopefile2.y scopefile2.z]-[scopefile1.x scopefile1.y scopefile1.z]);
stgshift = 1000*([scopefile2.x_mm scopefile2.y_mm scopefile2.z_mm]-[scopefile1.x_mm scopefile1.y_mm scopefile1.z_mm]);
pixshift = round(stgshift.*(dims-1)./imsize_um);
%%
% idaj : 1=right(+x), 2=bottom(+y), 3=below(+z)
% pixshift(iadj) = pixshift(iadj)+expensionshift(iadj); % initialize with a relative shift to improve CDP
%%
clc
projectionThr = 5;
debug = 0;
matchparams = modelParams(projectionThr,debug);
%%
iadj =find(gridshift);
if length(iadj)~=1 | max(iadj)>3
    error('not 6 direction neighbor')
end
[X_,Y_,out,rate_,pixshiftout,nonuniformity] = searchpair(desc1(:,1:3),desc2(:,1:3),pixshift,iadj,dims,matchparams);
if ~isempty(X_)
    X_ = correctTiles(X_,dims);
    Y_ = correctTiles(Y_,dims);
else
    matchparams_ = matchparams;
    matchparams_.opt.outliers = .5;
    [X_,Y_,out,rate_,pixshiftout,nonuniformity] = searchpair(desc1(:,1:3),desc2(:,1:3),pixshift,iadj,dims,matchparams_);
    if ~isempty(X_)
        X_ = correctTiles(X_,dims);
        Y_ = correctTiles(Y_,dims);
    end
end
%         paireddescriptor{ineig}.neigs = neigs_;
paireddescriptor.matchrate = rate_;
paireddescriptor.X = X_;
paireddescriptor.Y = Y_;
paireddescriptor.uni = mean(nonuniformity)<=.5;
%%
tag = 'XYZ';
if isempty(outfold)
    varargout{2} = paireddescriptor;
else
    outputfile = fullfile(outfold,sprintf('match-%s.mat',tag(iadj)));
    save(outputfile,'paireddescriptor','scopefile1','scopefile2')
    unix(sprintf('chmod g+rxw %s',outputfile))
end
end
function deployment(brain,tag,runlocal)
% mcc -m -v -R -singleCompThread /groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/zmatch_pipe/pointmatch.m -d /groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/pointmatch -a /groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/zmatch_pipe/thirdparty/CPD2 /groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/zmatch_pipe/functions
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
scopeparams = fullfile(matfolder,'scopeparams');
scopeparams = fullfile(matfolder,'scopeparams_pertile');

outfolder = fullfile(matfolder,'pointmatches')
mkdir(outfolder)

numcores = 16;
mkdir(fullfile(pwd,'shfiles'))
myfile = fullfile(pwd,'shfiles',sprintf('zmatchrun_%s_%s.sh',brain,date))
compiledfunc = '/groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/zmatch/zmatch'

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
% inds = 729:855;%round(linspace(729,855,127));

if 1
    myfiles = dir(fullfile(matfolder,'pointmatches','*.mat'));
    doneinds = cellfun(@(x) str2num(x(1:5)),{myfiles.name});
    [finished,bb] = min(pdist2((inds+1)',doneinds(:)),[],2);finished = ~finished;
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
    args = sprintf('''%s %s %s %s %s %05d %05d %0.1f %d %d''',compiledfunc,descriptorfile,scopefile,scopeparams,outfolder,inds(ii)+1,inds(ii+1),thr,numcores,zshift);
    %     mysub = sprintf('qsub -pe batch %d -l d_rt=%d -N %s -j y -o /dev/null -b y -cwd -V %s\n',numcores,esttime,name,args);
    mysub = sprintf('bsub -n%d -We %d -J %s -o /dev/null %s\n',numcores,esttime/60,name,args);
    if runlocal
        zmatch(sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/%s/matfiles/descriptors_ch0',[brain,tag]),...
            sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/%s/matfiles/scopeloc',[brain,tag]),...
            sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/%s/matfiles/scopeparams_pertile',[brain,tag]),...
            sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/%s/matfiles/pointmatches',[brain,tag]),...
            sprintf('%05d',inds(ii)+1), sprintf('%05d',inds(ii+1)), '0.1',num2str(feature('numCores')),'0')
    else
        fwrite(fid,mysub);
    end
end
unix(sprintf('chmod +x %s',myfile));
if ~runlocal
    fclose(fid);
end

end











