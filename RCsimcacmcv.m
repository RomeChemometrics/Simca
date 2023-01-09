function model=RCsimcacmcv(X,cl,modcl, nF, pX, cv,opt)

nin=find(cl==modcl);
Xin=X(nin,:);
nout=find(cl~=modcl);
Xout=X(nout,:);
clout=cl(nout);

cvtype=cv.cvtype;
segments=cv.cvsegments;

ntot=length(nin); %Total number of samples


if strcmp(cvtype, 'manual')
    manseg=segments;
    segments=max(manseg);
end

if strcmp(cvtype, 'loo')  %Loo cross-validation corresponds to syst123 with N segments
    cvtype='syst123';
    segments=ntot;
end

if strcmp(cvtype, 'random')  %Random cross-validation corresponds to syst123 with Random sorting
    rord=randperm(ntot);
end


acc_in=repmat({zeros(length(nin),1)}, 1,nF);
acc_out=repmat({zeros(length(nout),segments)}, 1,nF);
t2_in=repmat({zeros(length(nin),1)}, 1,nF);
t2_out=repmat({zeros(length(nout),segments)}, 1,nF);
q_in=repmat({zeros(length(nin),1)}, 1,nF);
q_out=repmat({zeros(length(nout),segments)}, 1,nF);
d_in=repmat({zeros(length(nin),1)}, 1,nF);
d_out=repmat({zeros(length(nout),segments)}, 1,nF);




for j=1:segments %Beginning of the crossvalidation loop
    
    switch cvtype
        case 'syst123'    %venetian blind cross-validation
            t=j:segments:ntot; %Validation set
            m=1:ntot; m(t)=[]; %Calibration set
            
        case 'random'    %Random cross-validation
            t=j:segments:ntot;
            t=rord(t); %Validation set
            m=1:ntot; m(t)=[]; %Calibration set
            
        case 'manual'    %Manual cross-validation
            t=find(manseg==j); %Validation set
            m=1:ntot; m(t)=[]; %Calibration set
            
        case 'syst111'    %contiguous blocks
            ns=ceil(ntot./segments);  %number of samples in each group
            if j==segments
                t=(j-1)*ns+1:ntot; %Validation set
            else
                t=(j-1)*ns+1:j*ns; %Calibration set
            end
            m=[1:(j-1)*ns ns*j+1:ntot]; %camp. nel set di calib.
    end
    
    Xm=Xin(m,:);
    Xt=Xin(t,:);
    
    for nn=1:nF
        
        sm=RCsimcacmcalc(Xm, nn, pX, opt);
        st=RCsimcacmapply(Xt, sm, opt);
        sout=RCsimcacmapply(Xout, sm, opt);
        
        acc_in{nn}(t)=st.accepted;
        acc_out{nn}(:,j)=sout.accepted;
        t2_in{nn}(t)=st.t2red;
        t2_out{nn}(:,j)=sout.t2red;
        q_in{nn}(t)=st.qred;
        q_out{nn}(:,j)=sout.qred;
        d_in{nn}(t)=st.dred;
        d_out{nn}(:,j)=sout.dred;
        
        
        
        
    end
    
end

nc=find(~isnan(clout)==1);
if ~isempty(nc)
    ncl=max(max(clout(nc)), modcl);
    ucl=unique(clout(nc))';
    spec=NaN(ncl, nF);
    tspec=zeros(1,nF);
    eff=zeros(1,nF);
    
end
sens=zeros(1,nF);



for i=1:nF
    sens(i)=100*sum(acc_in{i})/length(nin);
    if ~isempty(nc)
        tspec(i)=100*length(find(acc_out{i}(nc,:)==0))/(length(nc)*segments);
        
        for j=ucl
            spec(j,i)=100*length(find(acc_out{i}(nc(clout(nc)==j),:)==0))/(length(nc(clout(nc)==j))*segments);
        end
        
        eff(i)=sqrt(sens(i)*tspec(i));
    else
        tspec=[];
        spec=[];
        eff=[];
        
        
    end
    
    
end


model.accepted_modcl=acc_in;
model.accepted_others=acc_out;
model.t2red_modcl=t2_in;
model.t2red_others=t2_out;
model.qred_modcl=q_in;
model.qred_others=q_out;
model.dred_modcl=d_in;
model.dred_others=d_out;
model.samples_modcl=nin;
model.samples_others=nout;
model.class_others=clout;
model.sensitivity=sens;
model.totspecificity=tspec;
model.specificity=spec;
model.efficiency=eff;

switch opt.PCsel
    case 'auto'
        
        switch opt.PCcrit
            case 'compl'
                [~, optPC]=max(eff);
                
            case 'rig'
                [~, optPC]=find(sens>opt.minsens, 1, 'last');
                if isempty(optPC)
                    optPC=1; 
                end
                
        end
        
    case 'manual'
        
        figure('units', 'normalized', 'position', [0.2 0.2 0.5 0.5])
        if ~isempty(nc)
            plot(1:nF, sens, 'r', 'linewidth', 2.5)
            hold on
            plot(1:nF, tspec, 'b', 'linewidth', 2.5)
            plot(1:nF, eff, '--k', 'linewidth', 2.5)
            axis tight
            set(gca, 'fontsize', 16, 'fontweight', 'bold', 'linewidth', 2)
            xlabel('Number of PCs')
            ylabel('Figures of Merit')
            legend({'Sensitivity', 'Specificity', 'Efficiency'}, 'fontsize', 14, 'fontweight', 'bold', 'location', 'best')
        else
            plot(1:nF, sens, 'r', 'linewidth', 2.5)
            axis tight
            set(gca, 'fontsize', 16, 'fontweight', 'bold', 'linewidth', 2)
            xlabel('Number of PCs')
            ylabel('Figures of Merit')
            legend({'Sensitivity', 'Specificity', 'Efficiency'}, 'fontsize', 14, 'fontweight', 'bold', 'location', 'best')
        end
        
        optPC=input('Select the optimal number of PCs:   ');
        
end

model.optPC=optPC;
if ~isempty(nc)
    
    model.optsensitivity=sens(optPC);
    model.opttotspecificity=tspec(optPC);
    model.optspecificity=spec(:,optPC)';
    model.optefficiency=eff(optPC);
end