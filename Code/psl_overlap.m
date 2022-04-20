function [Min,Max,N] = psl_overlap(input,T,ratio,outname,includeT)
% first binarize,then overlap
% if include non-zero positive value,set T to be empty
foldertmp = fileparts(outname);
if ~isfolder(foldertmp)
    mkdir(foldertmp)
end
if ~exist('includeT','var') || isempty(includeT)
    includeT = 0;
end
if ~exist('T','var') || isempty(T)
    T = 0;
end
[Data,Header] = y_Read(input{1});
Data(:) = 0;
Header.dt = [16,0];
N = length(input);
M = zeros(prod(Header.dim),N);
for i = 1:N
    data = y_Read(input{i});
    M(:,i) = data(:);
end
    Mtmp = M;
    if includeT
    M(Mtmp < T) = 0;
    M(Mtmp >= T) = 1;
    else
    M(Mtmp <= T) = 0;
    M(Mtmp > T) = 1;
    end
Mconj = sum(M,2);
Max = max(Mconj(:));
data0 = Data;
data0(:) =  Mconj;
Min = round(ratio*N);
Mthr = Mconj;
Mthr(Mconj < Min) = 0;
dataThr = Data;
dataThr(:) = Mthr;
outname_raw = strjoin(strsplit(outname,'.'),'_raw.');
y_Write(data0,Header,outname_raw)
y_Write(dataThr,Header,outname)
end