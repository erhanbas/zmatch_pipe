function [X_,Y_,rate_,pixshift,nonuniformity] = searchpair_edge(des_center,des_cadj_ori,pixshiftinit,iadj,tile_size,matchparams)
%SEACHPAIR Summary of this function goes here
%
% [OUTPUTARGS] = SEACHPAIR(INPUTARGS) Explain usage here
%
% Inputs:
%   descent: Num_descriptor-by-Num_feature numerical array of the tile
%   descadjori: Num_descriptor-by-Num_feature numerical array of the
%   neighboring tile
%   pixshiftinit: initially estimated pixel shift
%   iadj: index of adjacant tile that have relative position shift w.r.t.
%   the current tile (?)
%   dims: tile size
%   matchparams: 
%   
% Outputs:
%
% Examples:
%
% Provide sample usage code here
%
% See also: List related files here

% $Author: base $	$Date: 2016/11/03 16:02:56 $	$Revision: 0.1 $
% Copyright: HHMI 2016
if isfield(matchparams, 'max_num_desc')
    total_num_descriptor = matchparams.max_num_desc;
else
    total_num_descriptor = 5000;
end
pixshift = pixshiftinit;
% pixshift = zeros(1,size(descent,2));
% pixshift(1:length(pixshiftinit)) = pixshiftinit;
%search
flag_stop = false;
iter = 0;
R = zeros(1,50);
% nonuniformity = zeros(1,10);
% clear nonuniformity

[X_,Y_,neigs_,rate_] = deal([]);
while ~flag_stop & iter<50% run a search
    %%
    iter = iter + 1;
    des_cadj = bsxfun(@plus, des_cadj_ori, [pixshift, 0]);
    % Why +3 ? 
    nbound = [max(pixshift(iadj), min(des_cadj(:,iadj))), ...
        min(tile_size(iadj), max(des_center(:,iadj))) + 3];
    X = des_center(des_center(:,iadj) > nbound(1) & des_center(:,iadj) < nbound(2),:);
    Y = des_cadj(des_cadj(:,iadj) > nbound(1) & des_cadj(:,iadj) < nbound(2),:);
    % If too much descriptor, prefer the long one
    % Maybe better solution: chop skeletons into short pieces and sample uniformly in
    % the overlapping space. 
    if size(X, 1) > total_num_descriptor
        tmp_idx_list = fun_bin_data_to_idx_list(X(:,4));
        tmp_seg_length = cellfun(@numel, tmp_idx_list);
        [tmp_seg_length, tmp_seg_idx] = sort(tmp_seg_length, 'descend');
        tmp_cumsum_voxel = cumsum(tmp_seg_length);
        [~, cutoff_idx] = min(abs(tmp_cumsum_voxel - total_num_descriptor));
        X = X(cat(2, tmp_idx_list{tmp_seg_idx(1:cutoff_idx)}), 1:3);
    else
        X = X(:, 1:3);
    end
    if size(Y, 1) > total_num_descriptor
        tmp_idx_list = fun_bin_data_to_idx_list(Y(:,4));
        tmp_seg_length = cellfun(@numel, tmp_idx_list);
        [tmp_seg_length, tmp_seg_idx] = sort(tmp_seg_length, 'descend');
        tmp_cumsum_voxel = cumsum(tmp_seg_length);
        [~, cutoff_idx] = min(abs(tmp_cumsum_voxel - total_num_descriptor));
        Y = Y(cat(2, tmp_idx_list{tmp_seg_idx(1:cutoff_idx)}), 1:3);
    else
        Y = Y(:, 1:3);
    end
    
    %%
    if size(X,1)<3 || size(Y,1)<3% not enough sample to match
        [X_,Y_,rate_,pixshift,nonuniformity] = deal([]);
        flag_stop = 1;
    else
        %% check uniformity of data
        nbins = [2 2];
        edges = {};
        for ii = 2:-1:1%length(dims)%[1 2 3],
            minx = 0;
            maxx = tile_size(ii);
            binwidth = (maxx - minx) / nbins(ii);
            edges{ii} = minx + binwidth*(0:nbins(ii));
        end
        % Get the bivaritive histogram for the location of the feature
        % points. These point are counted in four quarters of the x-y plane
        [accArr] = hist3([X(:,1:2);Y(:,1:2)],'Edges',edges);
        accArr = accArr(1:2,1:2);
        % Check if the four quarter are balanced by checking if the
        % diagonal quarters are larger than the average while the
        % anti-diagonal quaters are less than the average, or vice versa.
        if ~all(sum( accArr>mean(accArr(:)) ) & sum(accArr>mean(accArr(:)),2)' )  
            % non uniform over quad-representation
            nonuniformity(iter) = 1;
        else
            nonuniformity(iter) = 0;
        end
        
        try
            %% Match the descriptor
            [rate,X_,Y_,tY_] = descriptorMatchforz(X,Y,pixshift,iadj,matchparams);
            if size(X_,1)<3
                rate = 0; % overparametrized system
            end
            R(iter) = rate;
            if iter>1 & R(iter)-R(iter-1)<0
                % if the current matching is worse than the previous
                % matching, go back to the previous matching and stop the
                % iteration
                flag_stop = 1;
                X_ = X_t_1;
                Y_ = Y_t_1;
                rate = R(iter-1);
            else
                X_t_1 = X_;
                Y_t_1 = Y_;
                if rate<.95 && iadj ==3% no match
                    pixshift = pixshift + [0 0 5]; % expand more
                    flag_stop = 0;
                    error('increase shift')
                else % match found
                    flag_stop = 1;
                end
            end
            % store pairs
            rate_ = rate;
        catch
            X_ = [];
            Y_ = [];
            disp('error')
        end
    end
end
end
function [rate,X_,Y_,tY_] = descriptorMatchforz_relaxed(X,Y,pixshift,iadj,params)
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
% Y_(:,iadj) = Y_(:,iadj)- pixshift(iadj);% move it back to original location after CDP
Y_ = Y_- ones(size(Y_,1),1)*pixshift;% move it back to original location after CDP

% %%
% % displacement field between follows a field curve on x&y due to
% % optics and deformation curve due to tissue and cut force on z
% dispvec = X_-Y_;
% x = dispvec(:,iadj);
% if iadj==1 % x-neighbor
%     % x : x-displacement
%     % y : y-location
%     y = X_(:,2);
%     bw = [2 220];
% elseif iadj==2 % y-neighbor
%     % x : y-displacement
%     % y : x-location
%     y = X_(:,1);
%     bw=[3 100];
% else % z-neighbor
%     % x : z-displacement
%     % y : y-location (not too much on x as cut is on y direction)
%     y = X_(:,2);
%     bw=[2 220];
% end
% % build a probabilistic model of displacement vectors
% N = 101;
% gridx = linspace(min(x),max(x),N);
% gridy = linspace(min(y),max(y),N);
% [density,bw] = ksdensity2d([x y],gridx,gridy,bw);density=density'/max(density(:));
% [xmin,ix] = min(pdist2(x,gridx'),[],2);
% [ymin,iy] = min(pdist2(y,gridy'),[],2);
% idx = sub2ind([N,N],iy,ix);
% prob_inliers = density(idx)>max(density(idx))*.25;
% x_inline = x(prob_inliers,:);
% y_inline = y(prob_inliers,:);
% %
% % fit curve model
% [~,im] = max(density,[],2);
% sgn = -2*((max(im)==im(end) | max(im)==im(1))-.5);
% pinit = [median(y) sgn*1e-5 median(x)];
% warning off
% out = nlinfit(y_inline, x_inline, model, pinit,optimopts);
% warning on
% % outlier rejection based on parametric model
% xest = feval(model,out,y);
% outliers = abs(x-xest)>2;
% X_ = X_(~outliers,:);
% Y_ = Y_(~outliers,:);
% tY_ = tY_(~outliers,:);
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

function [rate,X_,Y_,tY_] = descriptorMatchforz(X,Y,pixshift,iadj,params)
%DESCRIPTORMATCH Summary of this function goes here
%
% [OUTPUTARGS] = DESCRIPTORMATCH(INPUTARGS) Explain usage here
%
% Inputs:
%   X, Y: two 2D real, double marices, specifying the position of the point
%   in the point cloud.
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
out = [];
opt = params.opt;
model = params.model;
optimopts = params.optimopts;
projectionThr = params.projectionThr;
debug = params.viz;
%% Initial match based on point drift
[Transform, C] = cpd_register(X,Y,opt);
%% check if match is found
% Compute the pairwise euclidean distance between two input array
pD = pdist2(X,Transform.Y);
[aa1,bb1] = min(pD,[],1);
[aa2,bb2] = min(pD,[],2);
keeptheseY = find([1:length(bb1)]'==bb2(bb1));
keeptheseX = bb1(keeptheseY)';

disttrim = aa1(keeptheseY)' < projectionThr;
X_ = X(keeptheseX(disttrim),:);
Y_ = Y(keeptheseY(disttrim),:);
tY_= Transform.Y(keeptheseY(disttrim),:);
% Rate is the ratio of the number of matched pair of distance less than
% projectionThr over the total number of matched pairs
rate = sum(disttrim)/length(disttrim);
% [pixshift rate]
% if rate < .5 % dont need to continue
%     [X_,Y_,out] = deal(0);
%     return
% end
%%
% Y_(:,iadj) = Y_(:,iadj)- pixshift(iadj);% move it back to original location after CDP
Y_ = bsxfun(@minus, Y_, pixshift);
% Y_ = Y_ - ones(size(Y_,1),1)*pixshift;% move it back to original location after CDP

% %%
% % displacement field between follows a field curve on x&y due to
% % optics and deformation curve due to tissue and cut force on z
% dispvec = X_-Y_;
% x = dispvec(:,iadj);
% if iadj==1 % x-neighbor
%     % x : x-displacement
%     % y : y-location
%     y = X_(:,2);
%     bw = [2 220];
% elseif iadj==2 % y-neighbor
%     % x : y-displacement
%     % y : x-location
%     y = X_(:,1);
%     bw=[3 100];
% else % z-neighbor
%     % x : z-displacement
%     % y : y-location (not too much on x as cut is on y direction)
%     y = X_(:,2);
%     bw=[2 220];
% end
% % build a probabilistic model of displacement vectors
% N = 101;
% gridx = linspace(min(x),max(x),N);
% gridy = linspace(min(y),max(y),N);
% [density,bw] = ksdensity2d([x y],gridx,gridy,bw);density=density'/max(density(:));
% [xmin,ix] = min(pdist2(x,gridx'),[],2);
% [ymin,iy] = min(pdist2(y,gridy'),[],2);
% idx = sub2ind([N,N],iy,ix);
% prob_inliers = density(idx)>max(density(idx))*.25;
% x_inline = x(prob_inliers,:);
% y_inline = y(prob_inliers,:);
% %
% % fit curve model
% [~,im] = max(density,[],2);
% sgn = -2*((max(im)==im(end) | max(im)==im(1))-.5);
% pinit = [median(y) sgn*1e-5 median(x)];
% warning off
% out = nlinfit(y_inline, x_inline, model, pinit,optimopts);
% warning on
% % outlier rejection based on parametric model
% xest = feval(model,out,y);
% outliers = abs(x-xest)>2;
% X_ = X_(~outliers,:);
% Y_ = Y_(~outliers,:);
% tY_ = tY_(~outliers,:);
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
%% Subfunctions
function [bin_cell_array, varargout] = fun_bin_data_to_idx_list(data)
% fun_bin_data_to_idx_list bin the data according to their values and
% output the corresponding index list
% Input: 
%   data: numerical vector
% Output: 
%   bin_cell_array: cell array, each cell constains a vector, whose
%   components are the indices of the component of data that have the same
%   value. 
%   varargout: unique data value 
% Author: Xiang Ji ( Department of Physics, UC San Diego )
% Nov 28, 2018

num_data = numel(data);
if ~issorted(data)
    [data, idx_list ]= sort(data, 'ascend');
else
    idx_list = 1 : num_data;
end

bin_size = 0;
bin_data = data(1);
bin_idx = zeros(1, round(num_data/2));
est_num_bin = 500;
bin_value_list = zeros(est_num_bin,1);
bin_value_list(1) = data(1);
bin_cell_array = cell(est_num_bin,1);
num_bin = 0;
for idx = 1 : num_data
    tmp_data = data(idx);
    
    if tmp_data == bin_data
        bin_size = bin_size + 1;
        bin_idx(bin_size) = idx_list(idx);
    else
        num_bin = num_bin + 1;
        bin_cell_array{num_bin} = bin_idx(1 : bin_size);
        bin_data = tmp_data;
        bin_value_list(num_bin + 1) = bin_data;
        bin_idx(1) = idx_list(idx);
        bin_size = 1;
    end
end
num_bin = num_bin + 1;
bin_cell_array{num_bin} = bin_idx(1 : bin_size);
bin_cell_array(num_bin + 1 : end) = [];
bin_value_list = bin_value_list(1 : num_bin);
if nargout > 1
    varargout{1} = bin_value_list;
end
end