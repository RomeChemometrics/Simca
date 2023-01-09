function model=RCsimcapred(varargin)
if nargin==3
    X=varargin{1};
    cl=varargin{2};
    simcamod=varargin{3};
elseif nargin==2
    X=varargin{1};
    cl=NaN(size(X,1),1);
    simcamod=varargin{2};
end

for i=1:length(simcamod.options.modclass)
    model.CModels{i}.Pred=RCsimcacpred(X,cl,simcamod.options.modclass(i), simcamod.CModels{i}.Model, simcamod.options);
    model.CModels{i}.modclass=simcamod.options.modclass(i);
end

model.options=simcamod.options; 
    