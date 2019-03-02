function pointmatch_task_local(brain,inputfolder,descriptorfolder,matchfolder,matfolder,directions,ch,runlocal)
%% deploys pointmatch function on array task

deployment(brain,inputfolder,descriptorfolder,matchfolder,matfolder,directions,ch,runlocal);
end

function deployment(brain,inputfolder,descriptorfolder,matchfolder,matfolder,directions,ch,runlocal);
%%
addpath(genpath('./functions'))
% addpath(genpath('./common'))
matlab_cluster_path = '/usr/local/matlab-2017a';
compiledfunc = '/groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/pointmatch_task/pointmatch_task';
[mypath,mycomp] = fileparts(compiledfunc);
taskwrapper = fullfile(pwd,'cluster2.sh');
taskscript = fullfile(mypath,sprintf('run_%s.sh',mycomp));

if ~exist(fileparts(compiledfunc),'dir');
    mkdir(fileparts(compiledfunc));
    mfilename_ = mfilename('fullpath');
    unix(sprintf('mcc -m -v -R -I %s -d %s -a %s',mfilename_,fileparts(compiledfunc),fullfile(fileparts(mfilename_),'functions')))
    %,fullfile(fileparts(mfilename_),'common')
end

scopefile = fullfile(matfolder,'scopeloc.mat');
if ~exist(matfolder,'dir')
    mkdir(matfolder)
end
if exist(scopefile,'file')
    load(scopefile,'scopeloc')
else
    newdash = 1; % set this to 1 for datasets acquired after 160404
    [scopeloc] = getScopeCoordinates(inputfolder,newdash);% parse from acqusition files
    save(fullfile(matfolder,'scopeloc.mat'),'neighbors','scopeloc','imsize_um','experimentfolder','inputfolder')
end
[neighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
directionMap = containers.Map({'-X','-Y','X','Y','-Z','Z'},[ 2, 3, 4, 5, 6, 7]);
%%
maxnumofdesc = 10e3;
checkversion=1;
if 0
    pixinit = zeros(size(neighbors,1),3);
    nummatches = zeros(size(neighbors,1),1);
else
    % load finished tile matches. find badly matched or missing tile pairs
    [regpts,featmap] = loadMatchedFeatures(scopeloc,matchfolder,directions,checkversion);
    % initalize missing tiles based on knn
    numthr = 50;
    [pixinit,nummatches] = initTiles(featmap,directions,scopeloc,numthr);
    % sum(nummatches)
end
badtiles = nummatches<numthr & ~isnan(neighbors(:,directionMap(directions)));

%%
% generate a txt file with all tile pairs to be matched with optional shift
% values
% (filename,from,to,numcores,exitcode)
outlistfile = fullfile(pwd,'shfiles',sprintf('outlistfile_%s_%s.txt',brain,date));
if ~runlocal;fid = fopen(outlistfile,'w');else;fid=0;end
try parfor_progress(0);catch;end
parfor_progress(length(badtiles))
parfor ii = 1:length(badtiles)
    if ~badtiles(ii)
        parfor_progress;
        continue
    end
    
    %%
    tile1 = fullfile(descriptorfolder,scopeloc.relativepaths{ii});
    acqusitionfolder1 = fileparts(scopeloc.filepath{ii});
    iineig = neighbors(ii,directionMap(directions));
    tile2 = fullfile(descriptorfolder,scopeloc.relativepaths{iineig});
    acqusitionfolder2 = fileparts(scopeloc.filepath{iineig});
    outfold =fullfile(matchfolder,scopeloc.relativepaths{ii});
    if exist(fullfile(outfold,sprintf('match-%s-1.mat','Z')),'file')
        parfor_progress;
        continue
    end
    if runlocal
        try
            pointmatch(tile1,tile2,acqusitionfolder1,acqusitionfolder2,outfold,pixinit(ii,:),ch,maxnumofdesc,0);
        catch
            sprintf('%d',ii)
        end
    else
        fprintf(fid,'%s %s %s %s %s %f %f %f %s %f\n',tile1,tile2,acqusitionfolder1,acqusitionfolder2,outfold,pixinit(ii,:),ch,maxnumofdesc);
    end
    parfor_progress;
end
parfor_progress(0);
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




