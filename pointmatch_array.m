function varargout = pointmatch_array(descriptorfolder,scopefile,from,to,pixshift,numcores,exitcode)
%% deploys pointmatch function on array task
if nargin<1
    brain = '2017-09-19';tag='';runlocal = 0;deployment(brain,tag,runlocal);return
end

tr=load(scopefile,'scopeloc');scopeloc=tr.scopeloc;
[neighbors] = buildNeighbor(scopeloc.gridix(:,1:3));
if ischar(pixshift)
    pixshift = eval(pixshift); % pass initialization
end

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

function deployment(brain,tag,runlocal)
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
function finished = checkmissing(outputfolder,inds)
myfiles = dir(fullfile(outputfolder,'*.mat'));
doneinds = cellfun(@(x) str2num(x(1:5)),{myfiles.name});
[finished,bb] = min(pdist2((inds+1)',doneinds(:)),[],2);finished = ~finished;
end






