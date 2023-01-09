function model=RCsimcacmod(X,cl,modcl, nF, pX, opt)

nin=find(cl==modcl);
Xin=X(nin,:);

model=RCsimcacmcalc(Xin, nF, pX, opt);
model=RCsimcacmapply(X, model, opt);

sens=100*length(find(model.accepted(nin)==1))/length(nin);
model.sensitivity=sens;
nout=find(cl~=modcl&~isnan(cl)==1);
if ~isempty(nout)
    tspec=100*length(find(model.accepted(nout)==0))/length(nout);
    model.totspecificity=tspec;
    ncl=max(max(cl(~isnan(cl)==1)), modcl);
    spec=NaN(1,ncl);
    clout=cl(nout);
    ucl=unique(clout)';
    
    for i=ucl
        spec(i)=100*length(find(model.accepted(nout(clout==i))==0))/length(nout(clout==i));
    end
    model.specificity=spec;
    model.efficiency=sqrt(sens*tspec);
else
    model.totspecificity=[];
    model.specificity=[];
    model.efficiency=[];
end




