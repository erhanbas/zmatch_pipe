function varargout = pointmatch_vessel(tile1,tile2,acqusitionfolder1,acqusitionfolder2,outfold,pixshift,ch,maxnumofdesc,exitcode)
% pointmatch_vessel finds the corresponding position of the same object in
% the adjacent tiles by the following steps:
% 1. Use the microscope stage position information to estimate the relative
% displacement between two tiles. Load the descriptors from
% vesselDescriptor. 
% 2. If two tiles both contain large vessels near the boundary and the two
% tiles are adjacent in z direction, use masked intensity based FFT
% registration to estiamte the translation transformation between two
% tiles. 
% ( Need to decide when to use the displacement from masked FFT to update
% the pixshift for the following computation)
% 3. Apply Coherent Point Drift algorithm to find the correspondence
% between the the skeleton of the vessel from two tiles. Function used:
% searchpair_vessel
% 4. If both tiles contain large vessels near the boundary, use the vessel
% edge for point cloud registration. Function used: fun_searchpair_vessel_edges
% 5. Record and save all the registration result. 
%
% Modified from Erhan Bas's pointmatch by Xiang Ji (xiangji.ucsd@gmail.com)
% Date: Dec 14, 2018



%% Complied file setting
% compiledfunc = '/groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/pointmatch/pointmatch';
% if ~exist(fileparts(compiledfunc),'dir')
%     mkdir(fileparts(compiledfunc));
%     mfilename_ = mfilename('fullpath');
%     % unix(sprintf('mcc -m -v -R -singleCompThread %s -d %s -a %s',mfilename_,fileparts(compiledfunc),fullfile(fileparts(mfilename_),'functions')))
%     unix(sprintf('mcc -m -v %s -d %s -a %s',mfilename_,fileparts(compiledfunc),fullfile(fileparts(mfilename_),'functions')))
%     unix(sprintf('chmod g+x %s',compiledfunc))
%     %,fullfile(fileparts(mfilename_),'common')
% end
%% Path setting
if ~isdeployed
    addpath(genpath('./functions'))
end
%% Default input
if nargin<1
    rawfolder = '/groups/mousebrainmicro/mousebrainmicro/data/';
    classifierfolder = '/nrs/mouselight/pipeline_output/';
    sample = '2018-08-15-skeltest';
    tileid1 = '/2018-08-18/00/00152';
    tileid2 = '/2018-08-18/00/00442';
    tile1 = fullfile(classifierfolder,sample,'/stage_2_descriptor_output',tileid1);
    tile2 = fullfile(classifierfolder,sample,'/stage_2_descriptor_output',tileid2);
    acqusitionfolder1 = fullfile(rawfolder,sample,'Tiling',tileid1);
    acqusitionfolder2 = fullfile(rawfolder,sample,'Tiling',tileid2);
end

if nargin<5
    outfold = tile1;
    pixshift = '[0 0 0]';
    ch='1';
    maxnumofdesc=5e3;
    exitcode = 0;
elseif nargin < 6
    pixshift = '[0 0 0]';
    ch='1';
    maxnumofdesc=5e3;
    exitcode = 0;
elseif nargin <7
    ch='1';
    maxnumofdesc=5e3;
    exitcode = 0;
elseif nargin <8
    maxnumofdesc=5e3;
    exitcode = 0;
elseif nargin <9
    exitcode = 0;
end
%% Other preparations
if ischar(pixshift)
    pixshift = eval(pixshift); % pass initialization
end
if ischar(maxnumofdesc)
    maxnumofdesc=str2double(maxnumofdesc);
end
if ischar(exitcode)
    exitcode=str2double(exitcode);
end
%% Use the scope position to estimate the pixel shift
if length(ch)>1
    ch_desc={ch(1),ch(2)};
else
    ch_desc={ch};
end
varargout{1} = exitcode;
tile_size = [1024,1536,251];
empty_pixel_size_xyz = [40, 30, 0];
projectionThr = 5;
debug = 0;

tag = 'XYZ';
scopefile1 = readScopeFile(acqusitionfolder1);
scopefile2 = readScopeFile(acqusitionfolder2);
imsize_um = [scopefile1.x_size_um,scopefile1.y_size_um,scopefile1.z_size_um];
% estimate translation
gridshift = ([scopefile2.x scopefile2.y scopefile2.z]-[scopefile1.x scopefile1.y scopefile1.z]);
iadj =find(gridshift);
% Stage shift in micron
stgshift = 1000*([scopefile2.x_mm scopefile2.y_mm scopefile2.z_mm]-[scopefile1.x_mm scopefile1.y_mm scopefile1.z_mm]);
% Conver the stage shift to voxel distance
if all(pixshift==0)
    pixshift = round(stgshift.*(tile_size-1)./imsize_um);
end
paireddescriptor.stage_pixshift = pixshift;
%%
% check if input exists
% For the compatibility of the original pipeline
[~, ~, tmp_ext] = fileparts(tile1);
if isempty(tmp_ext) %~strcmp(tmp_ext, '.mat')
    tmp = dir(fullfile(tile1, '*tor.mat'));
    tile1 = fullfile(tmp.folder, tmp.name);
    tmp = dir(fullfile(tile2, '*tor.mat'));
    tile2 = fullfile(tmp.folder, tmp.name);
end

if ~isfile(tile1) || ~isfile(tile2)
    rate_ = 0;
    X_skel = [];
    Y_skel = [];
    uni = 0;
else
    descriptor_1 = load(tile1);
    descriptor_2 = load(tile2);
    if (descriptor_1.record.exist_blv || descriptor_2.record.exist_blv)
        % 1. Test if the large vessel mask is in the overlapping region -
        % store the information.
        % 2. Load the raw data, use 2d masked fft registration to estimate the
        % translational transformation.
        % 3. Use the estimation to initialize the edge registration
        
        % Load images
        tic
        tile_image_1 = deployedtiffread(descriptor_1.record.fp_image);
        tile_image_2 = deployedtiffread(descriptor_2.record.fp_image);
        % Flip tiles
        tile_image_1 = flip(flip(tile_image_1, 1), 2);
        tile_image_2 = flip(flip(tile_image_2, 1), 2);
%% Mask FFT registration on XY direction
% Doesn't work if
% 1. Register entire 2D sections directly
% 2. Register entire 3D overlapping region directly
% 3. Register 2D overlapping section
% Because:
% 1. Narrow boundary, no obvious feature
% 2. Image destortion due to homography is maximize on two side of
% the image.
% Conclusiong: Only use masked FFT for matching in Z direction. 
%         % Use tile 2 as moving image - Work for Z direction stitching first
% 
%         % The pixshift is [x_shift, y_shift, z_shift], while the bounding
%         % box etc are in [y, x, z]. Flip the pixshift for intensity
%         % registration here. 
%         descriptor_valid_bbox_mmxx = descriptor_1.record.valid_bbox_mmxx;
%         descriptor_valid_bbox_mmxx(5) = descriptor_valid_bbox_mmxx(5) - 10;
%         pixshift_yxz = pixshift([2,1,3]);
%         tile_size_yxz = tile_size([2,1,3]);
%         mask_seach_expansion = 20;
%         overlap_bbox_1_mmxx = [1, 1, 1, tile_size_yxz];
%         overlap_bbox_2_mmxx = [1, 1, 1, tile_size_yxz];
%         if gridshift(iadj) > 0
%             overlap_bbox_1_mmxx(1:3) = max(overlap_bbox_1_mmxx(1:3), pixshift_yxz) - mask_seach_expansion;
%             overlap_bbox_2_mmxx(4:6) = tile_size_yxz - pixshift_yxz + mask_seach_expansion;
%         elseif gridshift(iadj) < 0
%             overlap_bbox_2_mmxx(1:3) = tile_size_yxz + pixshift_yxz - mask_seach_expansion;
%             overlap_bbox_1_mmxx(4:6) = tile_size_yxz + pixshift_yxz + mask_seach_expansion;
%         end
%         overlap_bbox_1_mmxx(1:3) = max(overlap_bbox_1_mmxx(1:3), descriptor_valid_bbox_mmxx(1:3));
%         overlap_bbox_2_mmxx(1:3) = max(overlap_bbox_2_mmxx(1:3), descriptor_valid_bbox_mmxx(1:3));
%         overlap_bbox_1_mmxx(4:6) = min(overlap_bbox_1_mmxx(4:6), descriptor_valid_bbox_mmxx(4:6));
%         overlap_bbox_2_mmxx(4:6) = min(overlap_bbox_2_mmxx(4:6), descriptor_valid_bbox_mmxx(4:6));
%         test_sec = 60;
%         test_image_1 = tile_image_1(overlap_bbox_1_mmxx(1):overlap_bbox_1_mmxx(4), overlap_bbox_1_mmxx(2):overlap_bbox_1_mmxx(5), test_sec);
%         test_image_2 = tile_image_2(overlap_bbox_2_mmxx(1):overlap_bbox_2_mmxx(4), overlap_bbox_2_mmxx(2):overlap_bbox_2_mmxx(5), test_sec);
%         [tmp_translation, tmp_c, ~, ~] = MaskedTranslationRegistration(test_image_1, test_image_2, ...
%             test_image_1 > 1.3e4, test_image_2 > 1.3e4);
% %         
% %         % How to convert the translation back to the original coordinate? 
%         pixshift_xy = [overlap_bbox_1_mmxx(2) - overlap_bbox_2_mmxx(2) + tmp_translation(1),...
%             overlap_bbox_1_mmxx(1) - overlap_bbox_2_mmxx(1) + tmp_translation(2)];
%         figure;
%         subplot(1,4,1)
%         imshow(image_sec_1);
%         title('Section from tile 1');
%         subplot(1,4,2)
%         imshow(image_sec_2);
%         title('Section from tile 2');
%         subplot(1,4,3:4)
%         RA = imref2d(image_sec_size, [1, image_sec_size(2)], [1, image_sec_size(1)]);
%         RB = imref2d(image_sec_size, [pixshift_xy(1),image_sec_size(2) + pixshift_xy(1)], [pixshift_xy(2),image_sec_size(1) + pixshift_xy(2)]);
%         imshowpair(imadjust(image_sec_1), RA, imadjust(image_sec_2), RB, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', 'green-magenta')
%         title(sprintf('Overlap after translation (x,y) = (%d, %d) ', pixshift_xy(1), pixshift_xy(2)))
%% Mask FFT registration on Z direction         
% The following registration only works for the Z direction 
% Does this part need to be improved for more robust registration?
% For example, the image for registration can be part of the
% section, which can be cropped according to the position of the
% edge estimated from the edge subscripts. 
% Also, multiple sections of registration can be run for more
% reliable estimation of the translational displacement. 
        
% The pixshift is [x_shift, y_shift, z_shift], while the bounding
% box etc are in [y, x, z]. Flip the pixshift for intensity
% registration here.
        descriptor_valid_bbox_mmxx = descriptor_1.record.valid_bbox_mmxx;
        pixshift_yxz = pixshift([2,1,3]);
        empty_pixel_shift = [0,0,0];
        empty_pixel_shift(iadj) = empty_pixel_size_xyz(iadj);
        empty_pixel_shift_yxz = empty_pixel_shift([2,1,3]);
%         empty_pixel_size_yxz = [0, 0, 0];
        tile_size_yxz = tile_size([2,1,3]);
        mask_seach_expansion = 0;
        overlap_bbox_1_mmxx = [1, 1, 1, tile_size_yxz];
        overlap_bbox_2_mmxx = [1, 1, 1, tile_size_yxz];
        if gridshift(iadj) > 0
            overlap_bbox_1_mmxx(1:3) = max(overlap_bbox_1_mmxx(1:3), pixshift_yxz + empty_pixel_shift_yxz) - mask_seach_expansion;
            overlap_bbox_2_mmxx(4:6) = tile_size_yxz - (pixshift_yxz) + mask_seach_expansion;
        elseif gridshift(iadj) < 0
            overlap_bbox_2_mmxx(1:3) = max(overlap_bbox_2_mmxx(1:3),  - (pixshift_yxz - empty_pixel_shift_yxz)) - mask_seach_expansion;
            overlap_bbox_1_mmxx(4:6) = tile_size_yxz + ( pixshift_yxz )+ mask_seach_expansion;
        end
        overlap_bbox_1_mmxx(1:3) = max(overlap_bbox_1_mmxx(1:3), descriptor_valid_bbox_mmxx(1:3));
        overlap_bbox_2_mmxx(1:3) = max(overlap_bbox_2_mmxx(1:3), descriptor_valid_bbox_mmxx(1:3));
        overlap_bbox_1_mmxx(4:6) = min(overlap_bbox_1_mmxx(4:6), descriptor_valid_bbox_mmxx(4:6));
        overlap_bbox_2_mmxx(4:6) = min(overlap_bbox_2_mmxx(4:6), descriptor_valid_bbox_mmxx(4:6));
        overlap_bbox_1_mmll = overlap_bbox_1_mmxx;
        overlap_bbox_1_mmll(4:6) = overlap_bbox_1_mmxx(4:6) - overlap_bbox_1_mmxx(1:3) + 1;
        overlap_bbox_2_mmll = overlap_bbox_2_mmxx;
        overlap_bbox_2_mmll(4:6) = overlap_bbox_2_mmxx(4:6) - overlap_bbox_2_mmxx(1:3) + 1;
%         3D Masked FFT
        est_int_th = 1.5e4;
        test_image_1 = crop_bbox3(tile_image_1, overlap_bbox_1_mmll, 'default');
        test_image_2 = crop_bbox3(tile_image_2, overlap_bbox_2_mmll, 'default');
        [tmp_translation, tmp_max_xcorr, tmp_c] = MaskedTranslationRegistration(test_image_1, test_image_2, ...
            test_image_1 > est_int_th , test_image_2 > est_int_th, [20,20,10]);
        fft_pixshift_xyz = overlap_bbox_1_mmll([2,1,3]) - overlap_bbox_2_mmll([2,1,3]) + tmp_translation';
        paireddescriptor.exist_blv = true;
        paireddescriptor.pixshift_mask_fft = fft_pixshift_xyz;
        paireddescriptor.matchrate_mask_fft = tmp_max_xcorr;
        toc
        % Visualization 
        vis_pixshift_xyz = fft_pixshift_xyz;
        vis_translation = vis_pixshift_xyz - overlap_bbox_1_mmll([2,1,3]) + overlap_bbox_2_mmll([2,1,3]);
        [test_image_2_moved]= imtranslate(test_image_2, vis_translation);
%         vis_image_2 = test_image_2(:, :, vis_sec - tmp_translation(3)); % Be careful about the minus sign. The z coordinate is pointing downward here. 
%         vis_image_2_moved = imtranslate(vis_image_2, tmp_translation(1:2));
        vis_sec = 50;
        vis_image_1 = test_image_1(:, :, vis_sec);
        vis_image_2 = test_image_2(:, :, vis_sec - vis_translation(3));
        figure;
        subplot(1,4,1);
        imshow(vis_image_1);
        title('Section from tile 1');
        subplot(1,4,2)
        imshow(vis_image_2);
        title('Section from tile 2');
        subplot(1,4,3)
        imshow(test_image_2_moved(:, :, vis_sec));
        title('Translated section from tile 2');
        subplot(1,4,4)
        imshowpair(vis_image_1, test_image_2_moved(:, :, vis_sec));
        title(sprintf('Image overlap: pixel shift (%d, %d, %d)', vis_pixshift_xyz));
%         paireddescriptor.mask_fft_mask_ratio = sec_mask_ratio;
        figure;
        subplot(1,4,1)
        imshow(tile_image_1(:, :, vis_sec));
        title('Section from tile 1');
        subplot(1,4,2)
        imshow(tile_image_2(:, :, vis_sec));
        title('Section from tile 2');
        subplot(1,4,3:4)
        RA = imref2d(tile_size_yxz(1:2), [1, tile_size_yxz(2)], [1, tile_size_yxz(1)]);
        RB = imref2d(tile_size_yxz(1:2), [vis_pixshift_xyz(1),tile_size_yxz(2) + vis_pixshift_xyz(1)], [vis_pixshift_xyz(2),tile_size_yxz(1) + vis_pixshift_xyz(2)]);
        imshowpair(imadjust(tile_image_1(:, :, vis_sec)), RA, imadjust(tile_image_2(:, :, vis_sec)), RB, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', 'green-magenta')
        title(sprintf('Overlap after translation (x,y) = (%d, %d %d)', vis_pixshift_xyz));
        clear tile_image_1 tile_image_2 test_image_1 test_image_2 
    else
        paireddescriptor.exist_blv = false;
        paireddescriptor.pixshift_mask_fft = [];
        paireddescriptor.matchrate_mask_fft = [];
        paireddescriptor.mask_fft_mask_ratio = [];
    end
    % Need to decide when to use this
    if ~isempty(paireddescriptor.pixshift_mask_fft)
%         pixshift = paireddescriptor.pixshift_mask_fft;
    end
%% Skeleton point registration
    if isempty(descriptor_1.skl_sub) || isempty(descriptor_2.skl_sub)
        rate_ = 0;
        X_skel = [];
        Y_skel = [];
        uni = 0;
    else
        desc1_skel = cat(2, correctTiles(descriptor_1.skl_sub,tile_size), descriptor_1.skl_label);
        desc2_skel = cat(2, correctTiles(descriptor_2.skl_sub,tile_size), descriptor_2.skl_label);
        matchparams = modelParams(projectionThr,debug); % Setting parameters for the matching algorithm
        matchparams.max_num_desc = maxnumofdesc;
        matchparams.scan_z_shift_Q = true;
        if length(iadj)~=1 || max(iadj)>3
            error('not 6 direction neighbor')
        end
        %% MATCHING
        [X_skel,Y_skel,rate_, pixshift_skl,nonuniformity] = searchpair_vessel(desc1_skel,desc2_skel,pixshift,iadj,tile_size,matchparams);
        if isempty(X_skel)
            % I am not sure if this step is very useful or not, sicne CPD
            % can drift points by quite a lot. 
            disp('No points found. Relax the outliers tolorance and search for matched pairs again');
            matchparams_ = matchparams;
            matchparams_.opt.outliers = .5;
            matchparams_.selected_close_descriptor_pair_Q = true;
            matchparams_.scan_z_shift_Q = false;
            [X_skel,Y_skel,rate_,pixshift_skl, nonuniformity] = searchpair_vessel(desc1_skel, desc2_skel, pixshift, iadj, tile_size, matchparams_);
        end
        
        if ~isempty(X_skel)
            X_skel = correctTiles(X_skel,tile_size);
            Y_skel = correctTiles(Y_skel,tile_size);
        end
        uni = mean(nonuniformity)<=.5;
        paireddescriptor.pixshift_skl = pixshift_skl;
        paireddescriptor.matchrate = rate_;
        paireddescriptor.X = X_skel;
        paireddescriptor.Y = Y_skel;
        paireddescriptor.uni = uni;
        if rate_ > 0.95
            pixshift = pixshift_skl;
        end
    end
%% If both the descriptor contains boundary large vessels 

    if (descriptor_1.record.compute_edge && descriptor_2.record.compute_edge)
        desc1_edge = descriptor_1.edge_sub;
        desc2_edge = descriptor_2.edge_sub;
        desc1_edge = correctTiles(desc1_edge, tile_size);
        desc2_edge = correctTiles(desc2_edge, tile_size);
        [X_edge, Y_edge, rate_edge, pixshift_edge] = fun_searchpair_vessel_edges(desc1_edge, desc2_edge, pixshift);
        if ~isempty(X_edge)
            X_edge = correctTiles(X_edge, tile_size);
            Y_edge = correctTiles(Y_edge, tile_size);
        end
        paireddescriptor.exist_blv = true;
        paireddescriptor.matchrate_edge = rate_edge;
        paireddescriptor.X_edge = X_edge;
        paireddescriptor.Y_edge = Y_edge;
        paireddescriptor.pixshift_edge = pixshift_edge;        
    else
        paireddescriptor.matchrate_edge = [];
        paireddescriptor.X_edge = [];
        paireddescriptor.Y_edge = [];
        paireddescriptor.pixshift_edge = [];
    end
end

if nargin>4
    if ~isempty(X_skel)
        %x:R, y:G, z:B
        col = median(Y_skel - X_skel, 1)+128;
        col = max(min(col,255),0);
        outpng = zeros(105,89,3);
        outpng(:,:,1) = col(1);
        outpng(:,:,2) = col(2);
        outpng(:,:,3) = col(3);
        if exist(fullfile(outfold,'Thumbs.png'),'file')
            % -f to prevent prompts
            unix(sprintf('rm -f %s',fullfile(outfold,'Thumbs.png')));
        end
        imwrite(outpng,fullfile(outfold,'Thumbs.png'))
        unix(sprintf('chmod g+rw %s',fullfile(outfold,'Thumbs.png')));
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

end
%% Subfunctions
function [Iout] = deployedtiffread(fileName,slices)
%DEPLOYEDTIFFREAD Summary of this function goes here
% 
% [OUTPUTARGS] = DEPLOYEDTIFFREAD(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: base $	$Date: 2015/08/21 12:26:16 $	$Revision: 0.1 $
% Copyright: HHMI 2015
warning off
info = imfinfo(fileName, 'tif');
if nargin<2
    slices = 1:length(info);
end
wIm=info(1).Width;
hIm=info(1).Height;
numIm = numel(slices);
Iout  = zeros(hIm, wIm, numIm,'uint16');

for i=1:numIm
    Iout(:,:,i) = imread(fileName,'Index',slices(i),'Info',info);
end

end