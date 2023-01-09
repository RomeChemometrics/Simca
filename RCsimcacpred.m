function model=RCsimcacpred(X,cl,modcl,cmodel, opt)

model=RCsimcacmapply(X, cmodel, opt);


nin=find(cl==modcl&~isnan(cl)==1);
if ~isempty(nin)

sens=100*length(find(model.accepted(nin)==1))/length(nin);
model.sensitivity=sens;
else
    model.sensitivity=[];
end

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




