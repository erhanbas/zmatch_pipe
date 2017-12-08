function varargout = pointmatch(tile1,tile2,acqusitionfolder1,acqusitionfolder2,outfold,pixshift,ch,maxnumofdesc,exitcode)
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
end

if nargin<5
    outfold = tile1;
    pixshift = '[0 0 0]';
    ch='1';
    maxnumofdesc=1e3;
    exitcode = 0;
elseif nargin < 6
    pixshift = '[0 0 0]';
    ch='1';
    maxnumofdesc=1e3;
    exitcode = 0;
elseif nargin <7
    ch='1';
    maxnumofdesc=1e3;
    exitcode = 0;
elseif nargin <8
    maxnumofdesc=1e3;
    exitcode = 0;
elseif nargin <9
    exitcode = 0;
end
if ischar(pixshift)
    pixshift = eval(pixshift); % pass initialization
end

if length(ch)>1
    ch_desc={ch(1),ch(2)};
else
    ch_desc={ch};
end
varargout{1} = exitcode;
dims = [1024,1536,251];
projectionThr = 5;
debug = 0;

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
% read descs
desc1 = readDesc(tile1,ch_desc);
desc2 = readDesc(tile2,ch_desc);

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
    desc1 = truncateDesc(desc1,maxnumofdesc);
    desc2 = truncateDesc(desc2,maxnumofdesc);
    % idaj : 1=right(+x), 2=bottom(+y), 3=below(+z)
    % pixshift(iadj) = pixshift(iadj)+expensionshift(iadj); % initialize with a relative shift to improve CDP
    matchparams = modelParams(projectionThr,debug);
    if length(iadj)~=1 | max(iadj)>3
        error('not 6 direction neighbor')
    end
    %% MATCHING
    [X_,Y_,out,rate_,pixshiftout,nonuniformity] = searchpair(desc1(:,1:3),desc2(:,1:3),pixshift,iadj,dims,matchparams);
    if isempty(X_)
        matchparams_ = matchparams;
        matchparams_.opt.outliers = .5;
        [X_,Y_,out,rate_,pixshiftout,nonuniformity] = searchpair(desc1(:,1:3),desc2(:,1:3),pixshift,iadj,dims,matchparams_);
    end
    
    if ~isempty(X_)
        X_ = correctTiles(X_,dims);
        Y_ = correctTiles(Y_,dims);
    end
    uni = mean(nonuniformity)<=.5;
end
paireddescriptor.matchrate = rate_;
paireddescriptor.X = X_;
paireddescriptor.Y = Y_;
paireddescriptor.uni = uni;
if nargin>4
    if ~isempty(X_)
        %x:R, y:G, z:B
        col = median(Y_-X_,3)+128;
        col = max(min(col,255),0);
        outpng = zeros(105,89,3);
        outpng(:,:,1) = col(1);
        outpng(:,:,2) = col(2);
        outpng(:,:,3) = col(3);
        imwrite(outpng,fullfile(outfold,'Thumbs-match.png'))
    end
%%
if isempty(outfold)
    varargout{2} = paireddescriptor;
else
    % if isempty(rate_); val=0;elseif rate_<1;val=0;else;val=1;end
    outputfile = fullfile(outfold,sprintf('match-%s.mat',tag(iadj))); % append 1 if match found
    %check if file exist
    if exist(outputfile,'file')
        % if main match exists, crete a versioned one
        outputfile1 = fullfile(outfold,sprintf('match-%s-1.mat',tag(iadj))); % append 1 if match found
        save(outputfile1,'paireddescriptor','scopefile1','scopefile2')
        unix(sprintf('chmod g+rw %s',outputfile1));
    else
        save(outputfile,'paireddescriptor','scopefile1','scopefile2')
        unix(sprintf('chmod g+rw %s',outputfile));
    end
end
end










