function varargout = pointmatch(tile1,tile2,acqusitionfolder1,acqusitionfolder2,outfold,pixshift,exitcode)
%%
if ~isdeployed
    addpath(genpath('./functions'))
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
    pixshift = '[0 0 0]'
    exitcode = 0;
elseif nargin < 6
    pixshift = '[0 0 0]'
    exitcode = 0;
elseif nargin <7
    exitcode = 0;
end
if ischar(pixshift)
    pixshift = eval(pixshift); % pass initialization
end
varargout{1} = exitcode;
dims = [1024,1536,251];
projectionThr = 5;
debug = 0;

%%
%%
tag = 'XYZ';
scopefile1 = readScopeFile(acqusitionfolder1);
scopefile2 = readScopeFile(acqusitionfolder2);
imsize_um = [scopefile1.x_size_um,scopefile1.y_size_um,scopefile1.z_size_um];
% estimate translation
gridshift = ([scopefile2.x scopefile2.y scopefile2.z]-[scopefile1.x scopefile1.y scopefile1.z]);
iadj =find(gridshift);
stgshift = 1000*([scopefile2.x_mm scopefile2.y_mm scopefile2.z_mm]-[scopefile1.x_mm scopefile1.y_mm scopefile1.z_mm]);
if all(pixshift==0)
    pixshift = round(stgshift.*(dims-1)./imsize_um);
end
%%
% check if valid output exists
if exist(fullfile(outfold,sprintf('match-%s-1.mat',tag(iadj))),'file')
    return
end

%%
% read descs
desc1 = readDesc(tile1,{'0'});
desc2 = readDesc(tile2,{'0'});

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
    desc1 = truncateDesc(desc1);
    desc2 = truncateDesc(desc2);
    % idaj : 1=right(+x), 2=bottom(+y), 3=below(+z)
    % pixshift(iadj) = pixshift(iadj)+expensionshift(iadj); % initialize with a relative shift to improve CDP
    clc
    matchparams = modelParams(projectionThr,debug);
    %%
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
    uni = mean(nonuniformity)<=.5;
end

paireddescriptor.matchrate = rate_;
paireddescriptor.X = X_;
paireddescriptor.Y = Y_;
paireddescriptor.uni = uni;

%%
if isempty(outfold)
    varargout{2} = paireddescriptor;
else
    outputfile = fullfile(outfold,sprintf('match-%s-%d.mat',tag(iadj),rate_>0)); % append 1 if match found
    save(outputfile,'paireddescriptor','scopefile1','scopefile2')
    unix(sprintf('chmod g+rxw %s',outputfile))
end
end










