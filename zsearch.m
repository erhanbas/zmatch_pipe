function zsearch(descriptorfile,scopefile,scopeparams,outfolder,indstart,indend,thr,numcores,xshift,yshift,zshift)
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
    brain = '2017-19-19';
    tag='';
    runlocal = 1;
    deployment(brain,tag,runlocal)
end
if ~isdeployed
    addpath(genpath('./thirdparty'))
    addpath(genpath('./functions'))
end
indstart = str2double(indstart);
indend = str2double(indend);
thr = str2double(thr);
numcores = str2double(numcores);
xshift = str2double(xshift);
yshift = str2double(yshift);
zshift = str2double(zshift);
expensionshift = [xshift yshift zshift];
% addpath(genpath('./thirdparty'))
outfile = fullfile(outfolder,sprintf('%05d_%05d-pointmatch',indstart,indend));
% load descriptor file
load(descriptorfile,'descriptors')
% load scopelocation file
load(scopefile, 'scopeloc', 'imsize_um', 'experimentfolder', 'inputfolder')
% load scopeparams file
load(scopeparams,'scopeparams')
[neighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
checkthese = [1 4 5 7]; % 0 - below
neigs = neighbors(indstart:indend,checkthese);

paireddescriptor = pointmatch(descriptors,neigs,scopeloc,scopeparams,thr,numcores,expensionshift);

save(outfile,'paireddescriptor')
end
%% Sub-function - deployment
function deployment(brain,tag,runlocal)
% mcc -m -v -R -singleCompThread /groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/stitching/zsearch.m -d /groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/zsearch -a /groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/stitching/thirdparty/CPD2 /groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/stitching/functions
% mcc -m -v -R -I /groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/stitching/zsearch.m -d /groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/zsearch -a /groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/stitching/thirdparty/CPD2/ -a /groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/stitching/functions
%%
% brain = '2015-06-19';
% tag = '_backup'
if nargin<3
    runlocal = 0;
end
experimentfolder = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/%s%s/',brain,tag)
matfolder = fullfile(experimentfolder,'matfiles/');
descriptorfile = fullfile(matfolder,'descriptors_ch0');
scopefile = fullfile(matfolder,'scopeloc');
scopeparams = fullfile(matfolder,'scopeparams_pertile');
outfolder = fullfile(matfolder,'pointmatches')
outsecondfolder = fullfile(matfolder,'pointmatches_seconditer')
mkdir(outfolder)
mkdir(outsecondfolder)

numcores = 8;
myfile = sprintf('zmatchrun_%s_%s.sh',brain,date)
compiledfunc = '/groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/zsearch/zsearch'

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
% inds = 729:855;%round(linspace(729,855,127));

inds = round(linspace(0,Ntiles,Ntiles+1));
if 0
    inds = round(linspace(0,Ntiles,1000+1));
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
% check previous runs
if exist(fullfile(matfolder,'regpts.mat'),'file')
    regpts = load(fullfile(matfolder,'regpts'),'regpts');regpts=regpts.regpts;
    %%
    pixinit =nan(length(regpts),3); 
    for ii=1:length(regpts)
        pixinit(ii,:) = median(regpts{ii}.X-regpts{ii}.Y);
    end
    %%
    inliers = find(all(isfinite(pixinit),2));
    anchors = scopeloc.gridix(inliers,1:3);
    queries = scopeloc.gridix(:,1:3);
    IDX = knnsearch(anchors,queries,'K',1,'distance',@distfun);%W=[1 1 100000]
    pixinit = pixinit(inliers(IDX),:);
end

%%
if ~runlocal
    fid = fopen(myfile,'w');
end
% zsearch(descriptorfile,scopefile,scopeparams,outfolder,indstart,indend)
% for ii=1:length(inds)
thr = .1;
[xshift,yshift,zshift] = deal(0);

parfor ii=1:length(inds)-1
    %%
    if any(inliers==ii)%finished(ii)
        continue
    end
    ii
    %generate random string
    randString = s( ceil(rand(1,sLength)*numRands) );
    name = sprintf('zm_%05d-%s',ii,randString);
    args = sprintf('''%s %s %s %s %s %05d %05d %0.1f %d %d %d %d''',compiledfunc,descriptorfile,scopefile,scopeparams,outfolder,inds(ii)+1,inds(ii+1),thr,numcores,xshift,yshift,zshift);
    mysub = sprintf('qsub -pe batch %d -l d_rt=%d -N %s -j y -o /dev/null -b y -cwd -V %s\n',numcores,esttime,name,args);
    if runlocal
        zsearch(descriptorfile,scopefile,scopeparams,...
            sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/%s/matfiles/pointmatches_seconditer',[brain,tag]),...
            sprintf('%05d',inds(ii)+1), sprintf('%05d',inds(ii+1)), ...
            '0.1','12',num2str(pixinit(ii,1)),num2str(pixinit(ii,2)),num2str(pixinit(ii,3)))
    else
        fwrite(fid,mysub);
    end
end
if ~runlocal
    unix(sprintf('chmod +x %s',myfile));
end

end

function paireddescriptor=pointmatch(descriptors,neigs,scopeloc,scopeparams,thr,numcores,expensionshift)
debug = 0;
if length(scopeparams)>1
    imsize_um = scopeparams(1).imsize_um;
else
    imsize_um = scopeparams.imsize_um;
end
projectionThr = 5; % distance between target and projected point has to be less than this number
dims = scopeparams.dims; % in xyz
% slid = [[75 960];[0 dims(2)];[0 dims(3)]];
% expensionshift = [0 0 zshift]; % HEURISTICS:: tissue expends, so overlap is bigger between tiles
%%
model = @(p,y) p(3) - p(2).*((y-p(1)).^2); % FC model
optimopts = statset('nlinfit');
optimopts.RobustWgtFun = 'bisquare';
% opt.method='nonrigid_lowrank';
opt.method='nonrigid';
opt.beta=6;            % the width of Gaussian kernel (smoothness)
opt.lambda=16;          % regularization weight
opt.viz=0;              % show every iteration
opt.outliers=0.9;       % use 0.7 noise weight
opt.fgt=0;              % do not use FGT (default)
opt.normalize=1;        % normalize to unit variance and zero mean before registering (default)
opt.corresp=1;          % compute correspondence vector at the end of registration (not being estimated by default)
%     opt.max_it=100;         % max number of iterations
%     opt.tol=1e-10;          % tolerance
matchparams.model = model;
matchparams.optimopts = optimopts;
matchparams.opt = opt;
matchparams.projectionThr = projectionThr;

matchparams.debug = ~isdeployed & debug;
matchparams.viz = ~isdeployed & 0;
%%
% checkthese = [1 4 5 7]; % 0 - right - bottom - below
% indicies are 1 based,e.g. x = 1:dims(1), not 0:dims(1)-1
% xyz_umperpix = zeros(size(neigs,1),3);
paireddescriptor = cell(size(neigs,1),1);
maxnumofdesc=3e3;
%%
% delete(gcp('nocreate'))
% parpool(numcores)
%%
for ineig = 1:size(neigs,1)
    %% load descriptor pairs X (center) - Y (adjacent tile)
    neigs_ = neigs(ineig,:);
    idxcent = neigs(ineig,1);
    
    % load center descriptors
    if 0
        descent = double(descriptors{idxcent}(:,1:3));
    elseif maxnumofdesc<size(descriptors{idxcent},1) & 0
        [vals,indssorted]=sort(descriptors{idxcent}(:,5),'descend'); % sort based on raw intensity
        validinds = indssorted(1:min(maxnumofdesc,length(indssorted)));
        descent = double(descriptors{idxcent}(validinds,1:3));
    else
        descent = double(descriptors{idxcent}(descriptors{idxcent}(:,5)>thr,1:3));
        descent(descent(:,3)>249,:)=[];
    end
    
    if size(descent,1)<3
        paireddescriptor{ineig}.X = [];
        paireddescriptor{ineig}.Y = [];
        paireddescriptor{ineig}.neigs = neigs_;
        paireddescriptor{ineig}.matchrate = 0;
        continue
    end
    
    for iadj = size(neigs,2)-1 %1:x-overlap, 2:y-overlap, 3:z-overlap
        %%
        idxadj =  neigs(ineig,iadj+1);
        if isnan(idxadj)
            paireddescriptor{ineig}.X = [];
            paireddescriptor{ineig}.Y = [];
            paireddescriptor{ineig}.neigs = neigs_;
            paireddescriptor{ineig}.matchrate = 0;
            continue
        end
        % load below descriptors
        if 0
            descadj = double(descriptors{idxadj}(:,1:3));
        elseif maxnumofdesc<size(descriptors{idxadj},1) & 0
            [vals,indssorted]=sort(descriptors{idxadj}(:,5),'descend');
            validinds = indssorted(1:min(maxnumofdesc,length(indssorted)));
            descadj = double(descriptors{idxadj}(validinds,1:3));
        else
            descadj = double(descriptors{idxadj}(descriptors{idxadj}(:,5)>thr,1:3));
            descadj(descadj(:,3)>249,:)=[];
        end
        if size(descadj,1)<3
            paireddescriptor{ineig}.X = [];
            paireddescriptor{ineig}.Y = [];
            paireddescriptor{ineig}.neigs = neigs_;
            paireddescriptor{ineig}.matchrate = 0;
            continue
        end
        %% initialize search based on scope location
        stgshift = 1000*(scopeloc.loc(neigs(ineig,1+iadj),:)-scopeloc.loc(neigs(ineig,1),:));
        pixshift = round(stgshift.*(dims-1)./imsize_um);
        % idaj : 1=right(+x), 2=bottom(+y), 3=below(+z)
        pixshift(iadj) = pixshift(iadj)+expensionshift(iadj); % initialize with a relative shift to improve CDP
        pixshift = expensionshift-[0 0 10];
        [X_,Y_,out,rate_,pixshiftout,nonuniformity] = pairsearch(descent,descadj,pixshift,iadj,dims,matchparams);
        %%
        initpixshift = -[[-41 114 -181];[29 120 -167];[-45 13 -186];[33 6 -174]];
        initgridloc = [[208 119	1287];[198	119	1287];[210 124	1287];[198 124 1287]];
        [~,sortedvals] = sort(pdist2(scopeloc.gridix(neigs_(1),1:2),initgridloc(:,1:2)),'ascend');
        initpixshift = initpixshift(sortedvals,:);
        
        for idx = 1: size(initpixshift,1)
            if isempty(X_)
                pixshift = initpixshift(idx,:);
                [X_,Y_,out,rate_,pixshiftout,nonuniformity] = pairsearch(descent,descadj,pixshift,iadj,dims,matchparams);
            end
        end
        
%         if isempty(X_)
%             pixshift = [378 652 215]-[338 786 37];
%             [X_,Y_,out,rate_,pixshiftout,nonuniformity] = pairsearch(descent,descadj,pixshift,iadj,dims,matchparams);
%         elseif isempty(X_)
%             pixshift = [441 957 251]-[392 1090 50];
%             [X_,Y_,out,rate_,pixshiftout,nonuniformity] = pairsearch(descent,descadj,pixshift,iadj,dims,matchparams);
%         end
        if isempty(X_) & isfinite(nonuniformity) & nonuniformity
            matchparams_ = matchparams;     
            matchparams_.opt.outliers = .5;
            [X_,Y_,out,rate_,pixshiftout,nonuniformity] = pairsearch(descent,descadj,pixshift,iadj,dims,matchparams_);
        end
        %mout(iadj,:) = out;
        paireddescriptor{ineig}.neigs = neigs_;
        paireddescriptor{ineig}.matchrate = rate_;
        paireddescriptor{ineig}.X = X_;
        paireddescriptor{ineig}.Y = Y_;
        paireddescriptor{ineig}.uni = mean(nonuniformity)<=.5;
    end
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

function [X_,Y_,mout_,rate_,pixshift,nonuniformity_] = ...
    pairsearch(descent,descadjori,pixshiftinit,iadj,dims,matchparams)

pixshift = pixshiftinit;
%search
flag = 0;
iter = 0;

clear res
nonuniformity = nan;
% check uniformity of data
nbins = [2 2];
edges = cell(1,2);
for ii=1:2%length(dims)%[1 2 3],
    minx = 0;
    maxx = dims(ii);
    binwidth = (maxx - minx) / nbins(ii);
    edges{ii} = minx + binwidth*(0:nbins(ii));
end

% [X_,Y_,mout_,neigs_,rate_] = deal([]);
while ~flag & iter<250 & all(pixshift<dims)% run a search
    iter
    nonuniformity(iter+1) = nan;
    % set to initialization
    descadj = descadjori + ones(size(descadjori,1),1)*pixshift;
%         % crop a slab
%     nbound = [0 0];
%     nbound(1) = max(pixshift(iadj),min(descadj(:,iadj)));
%     nbound(2) = min(dims(iadj),max(descent(:,iadj)))+3;
%     zslabcent = descent(:,iadj)>nbound(1)&descent(:,iadj)<nbound(2);
%     zslabadj = descadj(:,iadj)>nbound(1)&descadj(:,iadj)<nbound(2);
%     zslabcent = zslabcent&descent(:,1)>0&descent(:,2)>0&descent(:,1)<dims(1)&descent(:,2)<dims(2);
%     zslabadj = zslabadj&descadj(:,1)>0&descadj(:,2)>0&descadj(:,1)<dims(1)&descadj(:,2)<dims(2);
%     X = descent(zslabcent,:);
%     Y = descadj(zslabadj,:);

    % crop a slab
    [nbound] = zeros(3,2);
    zslabcent = ones(size(descent,1),1);
    zslabadj = ones(size(descadj,1),1);
    for idim=1:3
        nbound(idim,1) = max(0,max([pixshift(idim),min(descadj(:,idim)),min(descent(:,idim))]))-5;
        nbound(idim,2) = min([dims(idim),max(descent(:,idim)),max(descadj(:,idim))])+5;
        zslabcent = zslabcent & descent(:,idim)>nbound(idim,1) & descent(:,idim)<nbound(idim,2);
        zslabadj = zslabadj & descadj(:,idim)>nbound(idim,1) & descadj(:,idim)<nbound(idim,2);
    end
    X = descent(zslabcent,:);
    Y = descadj(zslabadj,:);
    
    if size(X,1)<3 | size(Y,1)<3% not enough sample to match
        rate = 0;
        X_ = [];
        Y_ = [];
        tY_ = [];
        out = zeros(1,3);
    else
        % check uniformity of data
        [accArr] = hist3([X(:,1:2);Y(:,1:2)],'Edges',edges);
        accArr = accArr(1:2,1:2);
        if ~all(sum(accArr>mean(accArr(:))) & sum(accArr>mean(accArr(:)),2)') % non uniform over quad-representation
            nonuniformity(iter+1) = 1;
        else
            nonuniformity(iter+1) = 0;
        end
        
        try
            %%
            [rate,X_,Y_,out,tY_] = descriptorMatchforz(X,Y,pixshift,iadj,matchparams);
            if 0
                figure(33)
                cla
                myplot3(X_,'.')
                myplot3(Y_,'.')
                myplot3(tY_,'ro')
            end
        catch
            rate = 0;
            X_ = [];
            Y_ = [];
            tY_ = [];
            out = zeros(1,3);
        end
        if size(X_,1)<3
            rate = 0; % overparametrized system
            X_ = [];
            Y_ = [];
            tY_ = [];
            out = zeros(1,3);
        end
    end
    res(iter+1).rate = rate;
    res(iter+1).X_ = X_;
    res(iter+1).Y_ = Y_;
    res(iter+1).out = out;
    res(iter+1).tY_ = tY_;
    res(iter+1).pixshift = pixshift;
    pixshift = pixshift + [0 0 5]; % expand more
    if rate>.95
        break
    end
    iter = iter + 1;
end
[rate,indmaxR] = max([res.rate]);
X_ = res(indmaxR).X_;
Y_ = res(indmaxR).Y_;
mout_ = res(indmaxR).out;
rate_ = res(indmaxR).rate;
pixshift = res(indmaxR).pixshift;
nonuniformity_ = nonuniformity(indmaxR);
end

function [rate,X_,Y_,out,tY_] = descriptorMatchforz(X,Y,pixshift,iadj,params)
%DESCRIPTORMATCH Summary of this function goes here
%
% [OUTPUTARGS] = DESCRIPTORMATCH(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2016/09/23 14:09:29 $	$Revision: 0.1 $
% Copyright: HHMI 2016
tY_ = [];
opt = params.opt;
model = params.model;
optimopts = params.optimopts;
projectionThr = params.projectionThr;
debug = params.viz;
%%
% initial match based on point drift
[Transform, C]=cpd_register(X,Y,opt);
%% check if match is found
pD = pdist2(X,Transform.Y);
[aa1,bb1]=min(pD,[],1);
[aa2,bb2]=min(pD,[],2);
keeptheseY = find([1:length(bb1)]'==bb2(bb1));
keeptheseX = bb1(keeptheseY)';

disttrim = aa1(keeptheseY)'<projectionThr;
X_ = X(keeptheseX(disttrim),:);
Y_ = Y(keeptheseY(disttrim),:);
tY_= Transform.Y(keeptheseY(disttrim),:);
rate = sum(disttrim)/length(disttrim);
% [pixshift rate]
if rate < .5 % dont need to continue
    [X_,Y_,out] = deal(0);
    return
end
%%
% BUG fix!!, pixshift is in 3D not only in z
% Y_(:,iadj) = Y_(:,iadj)- pixshift(iadj);% move it back to original location after CDP
Y_ = Y_- ones(size(Y_,1),1)*pixshift;% move it back to original location after CDP
%%
% displacement field between follows a field curve on x&y due to
% optics and deformation curve due to tissue and cut force on z
dispvec = X_-Y_;
x = dispvec(:,iadj);
if iadj==1 % x-neighbor
    % x : x-displacement
    % y : y-location
    y = X_(:,2);
    bw = [2 220];
elseif iadj==2 % y-neighbor
    % x : y-displacement
    % y : x-location
    y = X_(:,1);
    bw=[3 100];
else % z-neighbor
    % x : z-displacement
    % y : y-location (not too much on x as cut is on y direction)
    y = X_(:,2);
    bw=[2 220];
end
% build a probabilistic model of displacement vectors
N = 101;
gridx = linspace(min(x),max(x),N);
gridy = linspace(min(y),max(y),N);
[density,bw] = ksdensity2d([x y],gridx,gridy,bw);density=density'/max(density(:));
[xmin,ix] = min(pdist2(x,gridx'),[],2);
[ymin,iy] = min(pdist2(y,gridy'),[],2);
idx = sub2ind([N,N],iy,ix);
prob_inliers = density(idx)>max(density(idx))*.25;
x_inline = x(prob_inliers,:);
y_inline = y(prob_inliers,:);
%
% fit curve model
[~,im] = max(density,[],2);
sgn = -2*((max(im)==im(end) | max(im)==im(1))-.5);
pinit = [median(y) sgn*1e-5 median(x)];
out = nlinfit(y_inline, x_inline, model, pinit,optimopts);
% outlier rejection based on parametric model
xest = feval(model,out,y);
outliers = abs(x-xest)>2;
X_ = X_(~outliers,:);
Y_ = Y_(~outliers,:);
tY_ = tY_(~outliers,:);
if debug
    xgridest = feval(model,out,gridy);
    figure(100),
    subplot(2,2,iadj),cla
    imagesc(density,'Xdata',[gridx],'Ydata',[gridy])
    axis tight
    hold on,
    plot(x,y,'m.')
    plot(x_inline,y_inline,'mo')
    plot(x(~outliers),y(~outliers),'gd')
    plot(xgridest,gridy,'r-')
    
    subplot(2,2,4),
    cla
    hold on
    plot3(X(:,1),X(:,2),X(:,3),'b+')
    plot3(Y(:,1),Y(:,2),Y(:,3),'m.')
    plot3(Transform.Y(:,1),Transform.Y(:,2),Transform.Y(:,3),'ro')
    pause(.5)
    drawnow
    
end
end

