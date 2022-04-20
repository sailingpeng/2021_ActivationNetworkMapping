function psl_coactivation(seed,prefix,r,output_dir)
if ~exist('r','var') || isempty(r)
    r = 6;
end
if ~ischar(seed)
    seed = ['[[',num2str(aa(1)),', ',num2str(aa(2)),', ',num2str(aa(3)),']]'];
end
str1 = ['import sys; from neurosynth.analysis import network; ',...
    'sys.path.append(''/data/disk2/pengshaoling/code/python'');',...
    'from psl_dataset import Dataset; ',newline,...
    'dataset = Dataset.load(''/data/disk2/pengshaoling/Data/neurosynth-data-master/dataset_psl.pkl'');',...
    newline,'network.coactivation(dataset, seed = ''',seed,''', threshold=0.1, prefix=''',...
    prefix,''',r = ',num2str(r),',output_dir = ''',output_dir,''')',newline];
tmpfile = '/data/disk2/pengshaoling/downloads/pycode.py';
fid = fopen(tmpfile,'wt');
fprintf(fid,'%s\n',str1);
fclose(fid);
commandStr = ['python3.8 ',tmpfile];
system(commandStr);
delete(tmpfile)
end