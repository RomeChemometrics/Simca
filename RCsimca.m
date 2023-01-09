function model=RCsimca(varargin) 


if nargin==4
    X=varargin{1};
    cl=varargin{2};
    pX=varargin{3};
    opt=varargin{4};
elseif nargin==3
    X=varargin{1};
    cl=varargin{2};
    pX=varargin{3};
    
    opt.modclass='all'; 
    opt.maxPC=10; 
    opt.t2lim='Fdist';
    opt.t2cl=0.95;
    opt.qlim='jm';
    opt.qcl=0.95;
    opt.cmcrit='alt';
    opt.dcl=0.95;
    opt.cvtype='syst123'; 
    opt.cvsegments=5; 
    opt.PCcrit='compl'; 
    opt.PCsel='auto'; 
    opt.minsens=95; 
    
elseif nargin==1 &&strcmp(varargin{1}, 'options')
    
    opt.modclass='all'; 
    opt.maxPC=10; 
    opt.t2lim='Fdist';
    opt.t2cl=0.95;
    opt.qlim='jm';
    opt.qcl=0.95;
    opt.cmcrit='alt';
    opt.dcl=0.95;
    opt.cv.cvtype='syst123'; 
    opt.cv.cvsegments=5; 
    opt.PCcrit='compl'; 
    opt.PCsel='auto'; 
    opt.minsens=95; 
    model=opt;
    return
end


if ischar(opt.modclass) && strcmp(opt.modclass, 'all')
    opt.modclass=unique(cl(~isnan(cl)==1))'; 
end

if length(pX)==1
    pX=repmat(pX, 1,max(cl(~isnan(cl)==1))); 
end

if length(opt.cv)==1
    opt.cv=repmat(opt.cv, 1, max(cl(~isnan(cl)==1)));
end

if length(opt.maxPC)==1
    opt.maxPC=repmat(opt.maxPC, 1, max(cl(~isnan(cl)==1)));
end


for i=1:length(opt.modclass)
    if strcmp(opt.cv(opt.modclass(i)).cvtype, 'none')
        model.CModels{i}.CV=[];
        model.CModels{i}.Model=RCsimcacmod(X,cl,opt.modclass(i), opt.maxPC(opt.modclass(i)), pX{opt.modclass(i)}, opt);
        model.CModels{i}.modclass=opt.modclass(i);
        
    else
      model.CModels{i}.CV=RCsimcacmcv(X,cl,opt.modclass(i), opt.maxPC(opt.modclass(i)), pX{opt.modclass(i)},opt.cv(opt.modclass(i)),opt);
      model.CModels{i}.Model=RCsimcacmod(X,cl,opt.modclass(i), model.CModels{i}.CV.optPC, pX{opt.modclass(i)}, opt);
      model.CModels{i}.modclass=opt.modclass(i);
        
    end
end

model.options=opt; 



