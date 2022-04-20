function [emoTz,emoT] = psl_cs_fmri_conseed_seed2seed(dfold,cname,sfile,sfile2,pslcmd,...
    outputmask,kernel_parfor,dataset_info,outmat)
% 
% __________________________________________________________________________
% SUMMARY OF psl_cs_fmri_conseed_seed2seed
% 
% calculate FC between two seeds based on 1000 normative subjects from HCP or GSP.
% 
% SYNTAX:
% [emoTz,emoT] = psl_cs_fmri_conseed_seed2seed(dfold,cname,sfile,sfile2,pslcmd,...
%     outputmask,kernel_parfor,dataset_info,outmat)
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
%        (string, cell array of character vector)
%        the first set of seed images in .nii format
% 
% sfile2
%        (string, cell array of character vector)
%        the second set of seed images in .nii format
% 
% pslcmd
%        (string) 
%        fixed values, normally, don't change it 
%        'seed'
%        
% outputmask
%        (string) 
%        roi regions to to calculate seed-based functional connectivity, if "whole brain", set it as empty string '' 
%        
% kernel_parfor
%        (integer) 
%        number of kernels for parfor loop, it can be set at any value since it is not memory consumming.
%        
% dataset_info
%        (string) 
%        configuration information for 1000 normative connectomes.fixed values, normally, don't change it
%        dataset_info = 'dataset_info.mat'
% 
% outmat
%        (string) 
%        the name of the .mat file that store resultant variable "emoTz and emoT"
% __________________________________________________________________________
% OUTPUTS:
%     
% emoT: n*m matrix of r , n and m are the sizes of sfile and sfile2 respectively, e.g., emoT(i,j) is the r between seed i in sfile and seed j in sfile2
% emoTz: n*m matrix of fisher z , n and m are the sizes of sfile and sfile2 respectively
% __________________________________________________________________________
% COMMENTS:
% 
% detailed introduction of GSP dataset can be found in this article: "Brain
% Genomics Superstruct Project initial data release with structural,
% functional, and behavioral measures".this function was revised from
% function "cs_fmri_conseed_seed_tc" in LeadDBS toolbox (www.lead-dbs.org),
% which is far more sophisticated. add path of this toolbox (
% '/data/disk2/pengshaoling/code/toolbox/Lead_DBS') before run.
% 

tic

% if ~isdeployed
%     addpath(genpath('/autofs/cluster/nimlab/connectomes/software/lead_dbs'));
%     addpath('/autofs/cluster/nimlab/connectomes/software/spm12');
% end
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
    kernel_parfor = 6;
end
if ~exist('dataset_info','var')
    dataset_info = 'dataset_info.mat';
end
disp(['Connectome dataset: ',cname,'.']);
% ocname=cname;
if ~exist('outmat','var')
    outmat = '/data/disk2/pengshaoling/trial/trial_trial/outmat.mat';
end
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
        %         maskuseidx=ismember(dataset.vol.outidx,omaskidx);
    else
        omaskidx=dataset.vol.outidx;
        %         maskuseidx=1:length(dataset.vol.outidx);
    end
else
    omaskidx=dataset.vol.outidx; % use all.
    %     maskuseidx=1:length(dataset.vol.outidx);
end

% owasempty=0;
% if ~exist('outputfolder','var')
%     outputfolder=ea_getoutputfolder(sfile,ocname);
% %     owasempty=1;
% else
%     if isempty(outputfolder) % from shell wrapper.
%         outputfolder=ea_getoutputfolder(sfile,ocname);
% %         owasempty=1;
%     end
%     if ~strcmp(outputfolder(end),filesep)
%         outputfolder=[outputfolder,filesep];
%     end
% end

% if strcmp(sfile{1}(end-2:end),'.gz')
%     %gunzip(sfile)
%     %sfile=sfile(1:end-3);
%     usegzip=1;
% else
%     usegzip=0;
% end
for t = 1:size(sfile,1)
    seed{t,1}=ea_load_nii(ea_niigz(sfile{t,1}));
    sweights=seed{t,1}.img(dataset.vol.outidx);
    sweightmx=repmat(sweights,1,1);
    
    sweightidx{t,1}=find(sweights);
    sweightidxmx{t,1}=double(sweightmx(sweightidx{t,1},:));
    sweights(isnan(sweights))=0;
    sweights(isinf(sweights))=0; %
    
    sweights(abs(sweights)<0.0001)=0;
    sweights=double(sweights);
end
for s=1:size(sfile2,1)
    %         if size(sfile(1,:),2)>1
    %             dealingwithsurface=1;
    %         else
    %             dealingwithsurface=0;
    %         end
    %     for 1=1
    %         if exist(ea_niigz(sfile{t,1}),'file')
    seed2{s,1}=ea_load_nii(ea_niigz(sfile2{s,1}));
    
    %         else
    %             if size(sfile(1,:),2)==1
    %                 ea_error(['File ',ea_niigz(sfile{t,1}),' does not exist.']);
    %             end
    %             switch 1
    %                 case 1
    %                     sidec='l';
    %                 case 2
    %                     sidec='r';
    %             end
    %             seed{t,1}=dataset.surf.(sidec).space; % supply with empty space
    %             seed{t,1}.fname='';
    %             seed{t,1}.img(:)=0;
    %         end
    %         if ~isequal(seed{t,1}.mat,dataset.vol.space.mat)
    %             oseedfname=seed{t,1}.fname;
    %
    %             try
    %                 seed{t,1}=ea_conformseedtofmri(dataset,seed{t,1});
    %             catch
    %                 keyboard
    %             end
    %             seed{t,1}.fname=oseedfname; % restore original filename if even unneccessary at present.
    %         end
    
    %         [~,seedfn{s,1}]=fileparts(sfile{s,1});
    %         %         [~,seedfn2{s,1}]=fileparts(sfile2{s,1});
    %         if dealingwithsurface
    %             sweights=seed{t,1}.img(:);
    %         else
    sweights2=seed2{s,1}.img(dataset.vol.outidx);
    
    %             seed2_struc = cell2mat(seed2);
    %             Seed2_struc = {seed2_struc.img};
    %             sweights3 = cell2mat(cellfun(@(x) x(dataset.vol.outidx),Seed2_struc,'uniformoutput',false));
    %         end
    
    sweights2(isnan(sweights2))=0;
    sweights2(isinf(sweights2))=0; %
    
    sweights2(abs(sweights2)<0.0001)=0;
    sweights2=double(sweights2);
    
    %         sweights3(isnan(sweights3))=0;
    %         sweights3(isinf(sweights3))=0; %
    %
    %         sweights3(abs(sweights3)<0.0001)=0;
    %         sweights3=double(sweights3);
    
    %         try
    %             options=evalin('caller','options');
    %         end
    %         if exist('options','var')
    %             if strcmp(options.lcm.seeddef,'parcellation')
    %                 sweights=round(sweights);
    %             end
    %         end
    % assure sum of sweights is 1
    %sweights(logical(sweights))=sweights(logical(sweights))/abs(sum(sweights(logical(sweights))));
    
    sweightmx2=repmat(sweights2,1,1);
    
    sweightidx2{s,1}=find(sweights2);
    sweightidxmx2{s,1}=double(sweightmx2(sweightidx2{s,1},:));
    
    %         sweightmx3=repmat(sweights3,1,1);
    %         sweightidx3{s,1}=find(sweights3);
    %         sweightidxmx3{s,1}=double(sweightmx3(sweightidx3{s,1},:));
    
    % %     end
end
numseed=length(sfile);
% try
%     options=evalin('caller','options');
% end
% if exist('options','var')
%     if strcmp(options.lcm.seeddef,'parcellation') % expand seeds to define
%         ea_error('Command not supported for parcellation as input.');
%     end
% end
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
% numVoxUse = length(omaskidx);
% % init vars:
% for s=1:numseed
%     fX{s}=nan(numVoxUse,numSubUse);
%     rhfX{s}=nan(10242,numSubUse);
%     lhfX{s}=nan(10242,numSubUse);
% end
% fcr = zeros(numseed,1);
disp('Iterating through subjects...');
%     swidxl = {}; swidxr = {};
%     swidxmxl = {}; swidxmxr = {};
%     swidx = {};
%     swidxmx = {};
%     swidx2 = {};
%     swidxmx2 = {};
%     if size(sfile(i,:),2)>1
%         isROISurf = true;
%         swidxl = sweightidx{i,1}; swidxr = sweightidx{i,2};
%         swidxmxl = sweightidxmx{i,1}; swidxmxr = sweightidxmx{i,2};
%     else
%         isROISurf = false;
%         swidx2 = sweightidx2{i};
%         swidxmx2 = sweightidxmx2{i};
%     end
%     fXi = nan(numVoxUse,numSubUse);
if kernel_parfor > 1
    p = gcp('nocreate');
    if isempty(p)
        parpool(kernel_parfor);
    elseif p.NumWorkers < kernel_parfor
        delete(p)
        parpool(kernel_parfor)
    end
end
%     fXi2 = nan(1,numSubUse);
fXi3 = [];
for subj = 1:numSubUse % iterate across 1000 subjects, parfor
    
    mcfi = usesubjects(subj);
    disp(['Subject ', num2str(mcfi, '%04d'),'/',num2str(numSubUse,'%04d'),'...']);
    howmanyruns = ea_cs_dethowmanyruns(dataset,mcfi);
    %         thiscorr = nan(numVoxUse,howmanyruns);
    %         thiscorr2 = nan(1,howmanyruns);
    thiscorr3 = [];
    subIDVol = dataset.vol.subIDs{mcfi};
    for run=1:howmanyruns % load data
        gmtcstruc = load([dfoldvol,subIDVol{run+1}]);
        gmtc = single(gmtcstruc.gmtc);
        stc = [];
        for i=1:length(sweightidx) %iterate through each seed in sfile
            %     swidx = sweightidx{i};
            %     swidxmx = sweightidxmx{i};
            stc(i,:)=mean(gmtc(sweightidx{i},:).*repmat(sweightidxmx{i},1,size(gmtc,2)),1); % seed time course
            %             stc2=mean(gmtc(swidx2,:).*repmat(swidxmx2,1,size(gmtc,2)),1); % seed time course
        end
        %             swidx2 = sweightidx2{i};
        %             swidxmx2 = sweightidxmx2{i};
        stc3 = [];
        for ii = 1:length(sweightidx2)
            
            stc3(ii,:) = mean(gmtc(sweightidx2{ii},:).*repmat(sweightidxmx2{ii},1,size(gmtc,2)),1); % seed time course
        end
        %             thiscorr(:,run)=corr(stc',gmtc(maskuseidx,:)','type','Pearson');
        %             thiscorr2(:,run)=corr(stc',stc2','type','Pearson');
        thiscorr3(:,:,run)=corr(stc',stc3','type','Pearson');
    end
         thiscorr3z = atanh(thiscorr3);
   
    %         fXi(:,subj) = mean(thiscorr,2);
    %         fXi2(1,subj) = mean(thiscorr2);
    fXi3(:,:,subj) = mean(thiscorr3,3,'omitnan');
    fXi3z(:,:,subj) = mean(thiscorr3z,3,'omitnan');
end


%     fcr(numseed,1) = mean(fXi2);
emoT = mean(fXi3,3);
emoTz = mean(fXi3z,3);

% disp(['Seed ',num2str(i,'%03d'),'/',num2str(numseed,'%03d'),' done'])
if ~exist(outmat,'file')
    i = 1;
    save(outmat,'i')
end
save(outmat,'emoTz','emoT')
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
