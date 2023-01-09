function model=simcacompcrit(X, cl,pX,ntest,perct, opt)

fc={'perc', 'Fdist','Fdistrig','chi2' };
qc={'perc', 'jm','chi2box'};
dc={'sim', 'alt', 'ci'};
pc={'rig', 'compl'};

ncl=max(cl);

for nr=1:ntest
    cc=0;
    
    mm=[]; tt=[];
    for cln=1:ncl
        sc=find(cl==cln);
        nc=length(sc);
        nct=round(perct*nc);
        r=randperm(nc);
        tt=[tt; sc(r(1:nct))];
        mm=[mm; sc(r(nct+1:end))];
    end
    
    mm=sort(mm); tt=sort(tt);
    
    
    
    for i=1:4
        for j=1:3
            for k=1:3
                for l=1:2
                    opt.t2lim=fc{i};
                    opt.qlim=qc{j};
                    opt.cmcrit=dc{k};
                    opt.PCcrit=pc{l};
                    cc=cc+1;
                    
                    m=RCsimca(X(mm,:), cl(mm) , pX, opt);
                    p=RCsimcapred(X(tt,:), cl(tt) ,m);
                    for C=1:ncl
                        sens{C}(cc,nr)=p.CModels{C}.Pred.sensitivity;
                        tspec{C}(cc,nr)=p.CModels{C}.Pred.totspecificity;
                        spec{C}(cc,:, nr)=p.CModels{C}.Pred.specificity;
                        eff{C}(cc,nr)=p.CModels{C}.Pred.efficiency;
                        lv{C}(cc,nr)=p.CModels{C}.Pred.nPC;
                    end
                    
                    fcrit{cc,1}=fc{i};
                    qcrit{cc,1}=qc{j};
                    cmcrit{cc,1}=dc{k};
                    pccrit{cc,1}=pc{l};
                    
                end
            end
        end
    end
    
    opt.cmcrit='dd';
    opt.t2lim='chi2pom';
    opt.qlim='chi2pom';
    opt.PCcrit='rig';
    cc=cc+1;
    m=RCsimca(X(mm,:), cl(mm) , pX, opt);
    p=RCsimcapred(X(tt,:), cl(tt) ,m);
    for C=1:ncl
        sens{C}(cc, nr)=p.CModels{C}.Pred.sensitivity;
        tspec{C}(cc,nr)=p.CModels{C}.Pred.totspecificity;
        spec{C}(cc,:, nr)=p.CModels{C}.Pred.specificity;
        eff{C}(cc, nr)=p.CModels{C}.Pred.efficiency;
        lv{C}(cc,nr)=p.CModels{C}.Pred.nPC;
    end
    
    fcrit{cc,1}='chi2pom';
    qcrit{cc,1}='chi2pom';
    cmcrit{cc,1}='dd';
    pccrit{cc,1}='rig';
    
    opt.PCcrit='compl';
    
    cc=cc+1;
    m=RCsimca(X(mm,:), cl(mm) , pX, opt);
    p=RCsimcapred(X(tt,:), cl(tt) ,m);
    for C=1:ncl
        sens{C}(cc,nr)=p.CModels{C}.Pred.sensitivity;
        tspec{C}(cc,nr)=p.CModels{C}.Pred.totspecificity;
        spec{C}(cc,:, nr)=p.CModels{C}.Pred.specificity;
        eff{C}(cc,nr)=p.CModels{C}.Pred.efficiency;
        lv{C}(cc,nr)=p.CModels{C}.Pred.nPC;
    end
    
    fcrit{cc,1}='chi2pom';
    qcrit{cc,1}='chi2pom';
    cmcrit{cc,1}='dd';
    pccrit{cc,1}='compl';
    
end

for C=1:ncl
    model.CMres{C}.sensitivity=sens{C}; 
    model.CMres{C}.specificity=spec{C}; 
    model.CMres{C}.totspecificity=tspec{C};
    model.CMres{C}.efficiency=eff{C}; 
    model.CMres{C}.pc=lv{C}; 
    model.CMres{C}.meansensitivity=mean(sens{C},2); 
    model.CMres{C}.stdsensitivity=std(sens{C},[],2);
    model.CMres{C}.meanspecificity=mean(spec{C},3); 
    model.CMres{C}.stdspecificity=std(spec{C},[],3); 
    model.CMres{C}.meantotspecificity=mean(tspec{C},2);
    model.CMres{C}.stdtotspecificity=std(tspec{C}, [],2);
    model.CMres{C}.meanefficiency=mean(eff{C},2); 
    model.CMres{C}.stdefficiency=std(eff{C},[],2); 
    model.CMres{C}.meanpc=mean(lv{C},2); 
    model.CMres{C}.stdpc=std(lv{C},[],2); 
    model.CMres{C}.sensitivityres=[num2str(model.CMres{C}.meansensitivity, '%3.2f'), repmat(char(177), size(sens{C},1), 1), num2str(model.CMres{C}.stdsensitivity, '%3.2f') ]; 
    model.CMres{C}.totspecificityres=[num2str(model.CMres{C}.meantotspecificity, '%3.2f'), repmat(char(177), size(sens{C},1), 1), num2str(model.CMres{C}.stdtotspecificity, '%3.2f') ]; 
    
    for cc=1:ncl
        model.CMres{C}.specificityres{cc}=[num2str(model.CMres{C}.meanspecificity(:,cc), '%3.2f'), repmat(char(177), size(sens{C},1), 1), num2str(model.CMres{C}.stdspecificity(:,cc), '%3.2f') ];
    end
    
    model.CMres{C}.efficiencyres=[num2str(model.CMres{C}.meanefficiency, '%3.2f'), repmat(char(177), size(sens{C},1), 1), num2str(model.CMres{C}.stdefficiency, '%3.2f') ]; 
    model.CMres{C}.pcres=[num2str(model.CMres{C}.meanpc, '%3.2f'), repmat(char(177), size(sens{C},1), 1), num2str(model.CMres{C}.stdpc, '%3.2f') ]; 
    
    
    
end

    
    
model.combinations=[fcrit qcrit cmcrit pccrit]; 
model.options.ntest=ntest; 
model.options.perctest=perct; 

    
    
    
    
    
    
    
    
    
    