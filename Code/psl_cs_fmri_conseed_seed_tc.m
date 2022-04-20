function psl_cs_fmri_conseed_seed_tc(dfold,cname,sfile,pslcmd,...
    outputfolder,outputmask,kernel_parfor,dataset_info)
% 
% __________________________________________________________________________
% SUMMARY OF psl_cs_fmri_conseed_seed_tc
% 
% calculate seed based FC based on 1000 normative subjects from HCP or GSP.
% 
% SYNTAX:
%  psl_cs_fmri_conseed_seed_tc(dfold,cname,sfile,pslcmd,...
%        outputfolder,outputmask,kernel_parfor,dataset_info)
% __________________________________________________________________________
% INPUTS:
% 
% dfold
%        (string) 
%        the folder where 1000 normative connectomes were stored. fixed values, normally, don't change it 
%        '/data/disk2/pengshaoling/code/toolbox/Lead_DBS/connectomes'
%        
% cname
%        (string) 
%        two options: 'HCP1000>Full Set' or 'GSP1000>Full Set';
% 
% sfile
%        (string)
%        seed image in .nii format
% 
% pslcmd
%        (string) 
%        fixed values, normally, don't change it 
%        'seed'
%        
% outputfolder
%        (string) 
%        the folder to store results 
%        
% outputmask
%        (string) 
%        roi regions to to calculate seed-based functional connectivity, if "whole brain", set it as empty string '' 
%        
% kernel_parfor
%        (integer) 
%        number of kernels for parfor loop, for 1000 HCP dataset, the recommanded value is 1, otherwise it risks out of memory.
% 
%        
% dataset_info
%        (string) 
%        configuration information for 1000 normative connectomes.fixed values, normally, don't change it
%        dataset_info = 'dataset_info.mat'
% __________________________________________________________________________
% OUTPUTS:
%     
% ***_func_seed_AvgR.nii: mean r map;
% ***_func_seed_AvgR_Fz.nii: first fisher z, then mean
% ***_func_seed_T.nii: t map
% ***_func_seed_VarR.nii: variance map of r
% __________________________________________________________________________
% COMMENTS:
% 
% detailed introduction of GSP dataset can be found in this article: "Brain
% Genomics Superstruct Project initial data release with structural,
% functional, and behavioral measures".this function was revised from
% function "cs_fmri_conseed_seed_tc" from LeadDBS toolbox
% (www.lead-dbs.org). add path of this toolbox (
% '/data/disk2/pengshaoling/code/toolbox/Lead_DBS') before run. in
% xiaojiqun, better larger than 3, and in gaoxingneng better be 12,
% otherwise risk out of memory; for parfor in xiaojiqun, better 15 workers

tic
if ~iscell(sfile)
    sfile = {sfile};
end

if ~exist('dfold','var')
    dfold=''; % assume all data needed is stored here.
else
    if ~strcmp(dfold(end),filesep)
        dfold=[dfold,filesep];
    end
end
if ~exist('kernel_parfor','var')
    kernel_parfor = 1;
end
if ~exist('dataset_info','var')
    dataset_info = 'dataset_info.mat';
end
disp(['Connectome dataset: ',cname,'.']);
ocname=cname;

if ismember('>',cname)
    delim=strfind(cname,'>');
    subset=cname(delim+1:end);
    cname=cname(1:delim-1);
end

% prefs=ea_prefs;
% dfoldsurf=[dfold,'fMRI',filesep,cname,filesep,'surf',filesep];
dfoldvol=[dfold,'fMRI',filesep,cname,filesep,'vol',filesep]; % expand to /vol subdir.

d=load([dfold,'fMRI',filesep,cname,filesep,dataset_info]);
dataset=d.dataset;
clear d;
if exist('outputmask','var')
    if ~isempty(outputmask)
        omask=ea_load_nii(outputmask);
        omaskidx=find(omask.img(:));
        maskuseidx=ismember(dataset.vol.outidx,omaskidx);
    else
        omaskidx=dataset.vol.outidx;
        maskuseidx=1:length(dataset.vol.outidx);
    end
else
    omaskidx=dataset.vol.outidx; % use all.
    maskuseidx=1:length(dataset.vol.outidx);
end

if ~exist('outputfolder','var')
    outputfolder=ea_getoutputfolder(sfile,ocname);
else
    if isempty(outputfolder) % from shell wrapper.
        outputfolder=ea_getoutputfolder(sfile,ocname);
    end
    if ~strcmp(outputfolder(end),filesep)
        outputfolder=[outputfolder,filesep];
    end
end

if strcmp(sfile{1}(end-2:end),'.gz')
    %gunzip(sfile)
    %sfile=sfile(1:end-3);
    usegzip=1;
else
    usegzip=0;
end
for s=1:size(sfile,1)
    if size(sfile(s,:),2)>1
        dealingwithsurface=1;
    else
        dealingwithsurface=0;
    end
    if exist(ea_niigz(sfile{s,1}),'file')
        seed{s,1}=ea_load_nii(ea_niigz(sfile{s,1}));
    else
        if size(sfile(s,:),2)==1
            ea_error(['File ',ea_niigz(sfile{s,1}),' does not exist.']);
        end
        switch 1
            case 1
                sidec='l';
            case 2
                sidec='r';
        end
        seed{s,1}=dataset.surf.(sidec).space; % supply with empty space
        seed{s,1}.fname='';
        seed{s,1}.img(:)=0;
    end
    if ~isequal(seed{s,1}.mat,dataset.vol.space.mat) && (~dealingwithsurface)
        oseedfname=seed{s,1}.fname;
        
        try
            seed{s,1}=ea_conformseedtofmri(dataset,seed{s,1});
        catch
            keyboard
        end
        seed{s,1}.fname=oseedfname; % restore original filename if even unneccessary at present.
    end
    
    [~,seedfn{s,1}]=fileparts(sfile{s,1});
    if dealingwithsurface
        sweights=seed{s,1}.img(:);
    else
        sweights=seed{s,1}.img(dataset.vol.outidx);
    end
    sweights(isnan(sweights))=0;
    sweights(isinf(sweights))=0; %
    
    sweights(abs(sweights)<0.0001)=0;
    sweights=double(sweights);
    
    try
        options=evalin('caller','options');
    end
    if exist('options','var')
        if strcmp(options.lcm.seeddef,'parcellation')
            sweights=round(sweights);
        end
    end
    % assure sum of sweights is 1
    %sweights(logical(sweights))=sweights(logical(sweights))/abs(sum(sweights(logical(sweights))));
    sweightmx=repmat(sweights,1,1);
    
    sweightidx{s,1}=find(sweights);
    sweightidxmx{s,1}=double(sweightmx(sweightidx{s,1},:));
end

numseed=s;
try
    options=evalin('caller','options');
end
if exist('options','var')
    if strcmp(options.lcm.seeddef,'parcellation') % expand seeds to define
        ea_error('Command not supported for parcellation as input.');
    end
end

disp([num2str(numseed),' seeds, command = ',pslcmd,'.']);

numSubUse=length(dataset.vol.subIDs);

if ~exist('subset','var') % use all subjects
    usesubjects = 1:numSubUse;
else
    for ds=1:length(dataset.subsets)
        if strcmp(subset,dataset.subsets(ds).name)
            usesubjects=dataset.subsets(ds).subs;
            break
        end
    end
    numSubUse = length(usesubjects);
end

numVoxUse = length(omaskidx);

disp('Iterating through subjects...');
if kernel_parfor > 1
    p = gcp('nocreate');
    if isempty(p)
        parpool(kernel_parfor);
    elseif p.NumWorkers < kernel_parfor
        delete(p)
        parpool(kernel_parfor)
    end
end
fXi2 = single(zeros(length(sfile),numVoxUse,numSubUse));
parfor subj = 1:numSubUse % parfor iterate across 1000 normative subjects,
    mcfi = usesubjects(subj);
    disp(['Subject ', num2str(mcfi, '%04d'),'/',num2str(numSubUse,'%04d'),'...']);
    howmanyruns = ea_cs_dethowmanyruns(dataset,mcfi);
    %     thiscorr = nan(numVoxUse,howmanyruns);
    subIDVol = dataset.vol.subIDs{mcfi};
    
    thiscorr2 = [];
    for run=1:howmanyruns % load data for each run
        gmtcstruc = load([dfoldvol,subIDVol{run+1}]);
        gmtc = single(gmtcstruc.gmtc);
        
        %             stc=mean(gmtc(swidx,:).*repmat(swidxmx,1,size(gmtc,2)),1); % seed time course
        stc2 = [];% mean time series of voxels within seed.
        for ii = 1:length(sfile)
            stc2(ii,:)=mean(gmtc(sweightidx{ii},:).*repmat(sweightidxmx{ii},1,size(gmtc,2)),1);
        end
        %             thiscorr(:,run)=corr(stc',gmtc(maskuseidx,:)','type','Pearson');
        thiscorr2(:,:,run)=corr(stc2',gmtc(maskuseidx,:)','type','Pearson');
    end
    
    %         fXi(:,subj) = mean(thiscorr,2);
    fXi2(:,:,subj) = mean(thiscorr2,3);
end
parfor i = 1:length(sfile) % parfor
    fXi = squeeze(fXi2(i,:,:));
%     fXi2(1,:,:) = [];
%     if owasempty
%         outputfolder=ea_getoutputfolder(sfile(i),ocname);
%     end
    % export mean r map
    M=ea_nanmean(fXi',1);
    mmap=dataset.vol.space;
    mmap.dt=[16,0];
    mmap.img(:)=0;
    mmap.img=single(mmap.img);
    mmap.img(omaskidx)=M;
    
    mmap.fname=[outputfolder,seedfn{i},'_func_',pslcmd,'_AvgR.nii'];
    ea_write_nii(mmap);
    if usegzip
        gzip(mmap.fname);
        delete(mmap.fname);
    end
    
    % export variance map of r
    M=ea_nanvar(fXi');
    mmap=dataset.vol.space;
    mmap.dt=[16,0];
    mmap.img(:)=0;
    mmap.img=single(mmap.img);
    mmap.img(omaskidx)=M;
    
    mmap.fname=[outputfolder,seedfn{i},'_func_',pslcmd,'_VarR.nii'];
    ea_write_nii(mmap);
    if usegzip
        gzip(mmap.fname);
        delete(mmap.fname);
    end
    
    % fisher-transform:
    fXi=atanh(fXi);
    % export fz-mean, first fischer z, then average
    if any(isinf(fXi))
        warning(['seed ',sfile,' contains voxels with r = 1'])
    end
    M=nanmean(fXi');
    mmap=dataset.vol.space;
    mmap.dt=[16,0];
    mmap.img(:)=0;
    mmap.img=single(mmap.img);
    mmap.img(omaskidx)=M;
    mmap.fname=[outputfolder,seedfn{i},'_func_',pslcmd,'_AvgR_Fz.nii'];
    spm_write_vol(mmap,mmap.img);
    if usegzip
        gzip(mmap.fname);
        delete(mmap.fname);
    end
    % export T
    
    [~,~,~,tstat]=ttest(fXi');
    tmap=dataset.vol.space;
    tmap.img(:)=0;
    tmap.dt=[16,0];
    tmap.img=single(tmap.img);
    
    tmap.img(omaskidx)=tstat.tstat;
    
    tmap.fname=[outputfolder,seedfn{i},'_func_',pslcmd,'_T.nii'];
    spm_write_vol(tmap,tmap.img);
    if usegzip
        gzip(tmap.fname);
        delete(tmap.fname);
    end
    
    disp(['Seed ',num2str(i,'%03d'),'/',num2str(numseed,'%03d'),' done'])
end
toc


function s=ea_conformseedtofmri(dataset,s)
td=tempdir;
dataset.vol.space.fname=[td,'tmpspace.nii'];
ea_write_nii(dataset.vol.space);
s.fname=[td,'tmpseed.nii'];
ea_write_nii(s);

ea_conformspaceto([td,'tmpspace.nii'],[td,'tmpseed.nii']);
s=ea_load_nii(s.fname);
delete([td,'tmpspace.nii']);
delete([td,'tmpseed.nii']);


function howmanyruns=ea_cs_dethowmanyruns(dataset,mcfi)
if strcmp(dataset.type,'fMRI_matrix')
    howmanyruns=1;
else
    howmanyruns=length(dataset.vol.subIDs{mcfi})-1;
end
