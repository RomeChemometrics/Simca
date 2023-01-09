function model=RCsimcacmcalc(Xin, nF, pX, opt)


model=RCsimcapcamod(Xin, nF, pX);
ns=size(Xin, 1);


t2dof=[];
t2scfact=[];

qdof=[];
qscfact=[];

switch opt.t2lim
    case 'perc'
        [t2sort, ~]=sort(model.t2, 'ascend');
        qa=opt.t2cl*ns; ka=floor(qa);
        
        if ka~=0
            t2lim=t2sort(ka)+(qa-ka)*(t2sort(ka+1)-t2sort(ka));
        else
            t2lim=(qa-ka)*(t2sort(ka+1));
        end
        
    case 'Fdist'
        
        t2lim=nF*(ns-1)*finv(opt.t2cl,nF,ns-nF)./(ns-nF);
        t2dof=[nF,ns-nF];
        t2scfact=nF*(ns-1)/(ns-nF);
        
        
        
    case 'Fdistrig'
        t2lim=nF*((ns^2)-1)*finv(opt.t2cl,nF,ns-nF)./(ns*(ns-nF));
        t2dof=[nF,ns-nF];
        t2scfact=nF*((ns^2)-1)/(ns*(ns-nF));
        
    case 'chi2'
        t2lim=chi2inv(opt.t2cl, nF);
        t2dof=nF;
        t2scfact=1;
        
    case 'chi2pom'
        h0=mean(model.t2);
        Nh=max(round(2*(h0^2)/var(model.t2)), 1);
        t2lim=h0*chi2inv(opt.t2cl, Nh)/Nh;
        t2dof=Nh;
        t2scfact=h0;
end


switch opt.qlim
    
    case 'perc'
        [qsort, ~]=sort(model.q, 'ascend');
        qa=opt.qcl*ns; ka=floor(qa);
        
        if ka~=0
            qlim=qsort(ka)+(qa-ka)*(qsort(ka+1)-qsort(ka));
        else
            qlim=(qa-ka)*(qsort(ka+1));
        end
        
    case 'jm'
        
        theta1 = sum(model.detail.eigs(nF+1:end));
        theta2 = sum(model.detail.eigs(nF+1:end).^2);
        theta3 = sum(model.detail.eigs(nF+1:end).^3);
        if theta1==0
            qlim = 0;
        else
            h0= 1-((2*theta1*theta3)/(3*(theta2.^2)));
            if h0<0.001
                h0 = 0.001;
            end
            ca    = sqrt(2)*erfinv(2*opt.qcl-1);
            h1    = ca*sqrt(2*theta2*h0.^2)/theta1;
            h2    = theta2*h0*(h0-1)/(theta1.^2);
            qlim = theta1*(1+h1+h2).^(1/h0);
        end
        
    case 'chi2box'
        
        theta1 = sum(model.detail.eigs(nF+1:end));
        theta2 = sum(model.detail.eigs(nF+1:end).^2);
        
        
        g=theta2/theta1;
        Ng=(theta1^2)/theta2;
        qlim=g*chi2inv(opt.qcl, Ng);
        qdof=Ng;
        qscfact=g;
        
    case 'chi2pom'
        v0=mean(model.q);
        Nv=max(round(2*(v0^2)/var(model.q)), 1);
        qlim=v0*chi2inv(opt.qcl, Nv)/Nv;
        qdof=Nv;
        qscfact=v0;
        
end

model.t2lim=t2lim;
model.t2dof=t2dof;
model.t2scfact=t2scfact;

model.qlim=qlim;
model.qdof=qdof;
model.qscfact=qscfact;


switch opt.cmcrit
    
    case 'sim'
        dlim=1;
        
    case 'alt'
        
        dlim=sqrt(2);
        
    case 'ci'
        
        
        theta1 = sum(model.detail.eigs(model.nPC+1:end));
        theta2 = sum(model.detail.eigs(model.nPC+1:end).^2);
        tr1=(model.nPC/model.t2lim)+(theta1/model.qlim);
        tr2=(model.nPC/(model.t2lim^2))+(theta2/(model.qlim^2));
        gd=tr2/tr1;
        hd=(tr1^2)/tr2;
        dlim=gd*chi2inv(opt.dcl, hd);
        
    case 'dd'
        
        dlim=chi2inv(opt.dcl, model.t2dof+model.qdof);
        
end

model.dlim=dlim;

