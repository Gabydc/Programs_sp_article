function  savfig(dir,file,h,varargin)
% fprintf('Total number of inputs = %d\n',nargin);
nVarargs = length(varargin);
if nVarargs>0
 B=[dir   file varargin{1} '.fig'];
 saveas(h(1),B)
 B=[dir   file varargin{1} '.jpg'];
 saveas(h(1),B) 
else
   B=[dir   file  '.fig'];
 saveas(h(1),B)
 B=[dir   file '.jpg'];
 saveas(h(1),B) 
end
