function psl_vlsm(scan1,scan2,outdir,perm,mask,thr,if2tailed)
% outdir: no need to create outdir before this function was run
% psl_vlsm(F.tmap_vw,F.A.tmap_vw,F.f_vlsm) scan1 and scan2 must be
% binalized image "t_0.001_thresh.nii" is key image, by default, it is
% thresholded by corrected P = 0.05 with cluster-forming p(uncorrected) =
% 0.001; minimum 10% of subjects lesioned at a voxel for that voxel to be
% included in the analysis,
% "t_thresh.nii" is voxelwise FWE thresholded image

ratio = 0.1;
if ~exist('perm','var') || isempty(perm)
    perm = 5000;
end
if ~exist('if2tailed','var') || isempty(if2tailed)
    if2tailed = 1;
end
if if2tailed
    cl_p = 0.001/2;
else
    cl_p = 0.001;
end

formatOut = 'mmdd_HHMMSSFFF';
logdate = datestr(now,formatOut);

addpath(genpath([getenv('psldir'),'/code/toolbox/vlsm2']))
imgdir = [getenv('psldir'),'/downloads/',logdate,'tmp_vlsm_image'];
mkdir(imgdir)
cellfun(@(x) copyfile(x,imgdir),[scan1;scan2])
[~,D.txt1] = cellfun(@fileparts,scan2,'uniformoutput',false);
D.scan1 = strcat(imgdir,'/',D.txt1,'.nii');
[~,D.txt2] = cellfun(@fileparts,scan1,'uniformoutput',false);
D.scan2 = strcat(imgdir,'/',D.txt2,'.nii');
if exist('thr','var') && ~isempty(thr)
    cellfun(@(x,y) psl_binalizeimage(x,y,thr,1),[D.scan1;D.scan2],[D.scan1;D.scan2])
end

% create txt file
fileID = fopen([imgdir,'/design.txt'],'w');
fprintf(fileID,'%s\t%s\n','Filename','group');
Cell = [[D.txt1;D.txt2],num2cell([ones(length(D.txt1),1);zeros(length(D.txt2),1)])];
for i = 1:size(Cell,1)
    fprintf(fileID,'%s\t%d\n',Cell{i,:});
end
fclose(fileID);
if isfolder(outdir)
    rmdir(outdir,'s')
end
vlsm2([imgdir,'/design.txt'], imgdir, outdir,'nperms',perm,...
    'maskthresh',round(ratio*size(Cell,1)),'mask',mask,'alpha',cl_p)
rmpath(genpath([getenv('psldir'),'/code/toolbox/vlsm2']))
rmdir(imgdir,'s')