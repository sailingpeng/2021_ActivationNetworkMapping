function psl_SnPM_ttest1(scan,outdir,perm,mask,ifvw,if2tailed)
% perform one sample t text for specificity analysis in lesion network
% mapping, positive T means scan1 is larger than 0, default cluster
% defining threshold of P=0.001, cluster-wise FWE-corrected P=0.05). scan
% = g_ls([folder,filesep,'fmri',filesep,'*_PosT.nii']) see brain paper:
% "Network localization of heterogeneous neuroimaging findings"
% descriptioin of resultant files:
% http://www.nisox.org/Software/SnPM13/man,the key image is
% 'SnPM_filtered.nii', which is the t map, not z map.
rmpath(genpath([getenv('psldir'),'/code/toolbox/vlsm2']))
matlabbatch = {};
if ~exist('outdir','dir')
    mkdir(outdir)
end
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

[~,foldername] = fileparts(outdir);
if nargin == 2
    perm = 5000;
end
if ~exist('mask','var')
    mask = '/data/disk2/pengshaoling/fusion/Masks/Probability_2mm_tri_02.nii';
end
matlabbatch{1}.spm.tools.snpm.des.OneSampT.DesignName = ...
    'MultiSub: One Sample T test on diffs/contrasts';
matlabbatch{1}.spm.tools.snpm.des.OneSampT.DesignFile = 'snpm_bch_ui_OneSampT';
matlabbatch{1}.spm.tools.snpm.des.OneSampT.dir = {outdir};
matlabbatch{1}.spm.tools.snpm.des.OneSampT.P = strcat(scan,',1');
matlabbatch{1}.spm.tools.snpm.des.OneSampT.cov = struct('c', {}, 'cname', {});
matlabbatch{1}.spm.tools.snpm.des.OneSampT.nPerm = perm;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.vFWHM = [0 0 0];
matlabbatch{1}.spm.tools.snpm.des.OneSampT.bVolm = 1;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.ST.ST_U = cl_z;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.masking.tm.tm_none = 1;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.masking.im = 0;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.masking.em = {[mask,',1']};
matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalc.g_omit = 1;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalm.glonorm = 1;
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