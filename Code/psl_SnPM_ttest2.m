function psl_SnPM_ttest2(scan1,scan2,outdir,perm,mask,ifvw,if2tailed)
% perform two sample t text for specificity analysis in lesion network
% mapping, positive T means scan1 is larger than scan2, default cluster
% defining threshold of P=0.001, cluster-wise FWE-corrected P=0.05). scan1
% = g_ls([folder,filesep,'fmri',filesep,'*_PosT.nii']) see brain paper:
% "Network localization of heterogeneous neuroimaging findings"
% descriptioin of resultant files:
% http://www.nisox.org/Software/SnPM13/man,
%the key image is 'SnPM_filtered.nii', which is the t map, not z map. its
%corresponding unthresholded image is "snpmT+.nii"
% default is two-tailed
%-----------------------------------------------------------------------
rmpath(genpath([getenv('psldir'),'/code/toolbox/vlsm2']))
if ~exist('ifvw','var') || isempty(ifvw)
    ifvw = 0;
end
if ~exist('if2tailed','var') || isempty(if2tailed)
    if2tailed = 1;
end
if if2tailed
    cl_z = 3.29;
    vw_p = 0.05/2;
else
    cl_z = 3.09023230616781;
    vw_p = 0.05;
end
matlabbatch = {};
if ~exist('outdir','dir')
    mkdir(outdir)
end
[~,foldername] = fileparts(outdir);
if nargin == 2
    perm = 5000;
end
if ~exist('mask','var')
    mask = '/data/disk2/pengshaoling/fusion/Masks/Probability_2mm_tri_02.nii';
end
Cell = {strcat(scan1,',1'),strcat(scan2,',1')};
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.DesignName = ...
    '2 Groups: Two Sample T test; 1 scan per subject';
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.DesignFile = 'snpm_bch_ui_TwoSampT';
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.dir = {outdir};
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.scans1 = Cell{1};
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.scans2 = Cell{2};
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov = struct('c', {}, 'cname', {});
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.nPerm = perm;
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.vFWHM = [0 0 0];
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.bVolm = 1;
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.ST.ST_U = cl_z;
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.masking.tm.tm_none = 1;
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.masking.im = 0;
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.masking.em = {[mask,',1']};
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalc.g_omit = 1;
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalm.glonorm = 1;
matlabbatch{2}.spm.tools.snpm.cp.snpmcfg = {[outdir,filesep,'SnPMcfg.mat']};
matlabbatch{3}.spm.tools.snpm.inference.SnPMmat = {[outdir,filesep,foldername,'.mat']};
if ifvw
    matlabbatch{3}.spm.tools.snpm.inference.Thr.Vox.VoxSig.FWEth = vw_p;
else
    matlabbatch{3}.spm.tools.snpm.inference.Thr.Clus.ClusSize.CFth = NaN;
    matlabbatch{3}.spm.tools.snpm.inference.Thr.Clus.ClusSize.ClusSig.FWEthC = 0.05;
end
matlabbatch{3}.spm.tools.snpm.inference.Tsign = 1;
matlabbatch{3}.spm.tools.snpm.inference.WriteFiltImg.name = [foldername,'_filtered'];
matlabbatch{3}.spm.tools.snpm.inference.Report = 'MIPtable';
spm('defaults', 'FMRI');
spm_jobman('run',matlabbatch);
% clear matlabbatch
List = g_ls([outdir,filesep,'*.hdr']);
cellfun(@psl_hdr2nii,List);
Listimg = g_ls([outdir,filesep,'*.img']);
listdel = [List;Listimg];
cellfun(@delete,listdel)
end