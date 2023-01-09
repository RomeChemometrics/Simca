function model=RCsimcacmapply(Xout, cmod, opt)

ns=size(Xout,1);

model=RCsimcapcapred(Xout, cmod); 
model.t2lim=cmod.t2lim;
model.t2dof=cmod.t2dof; 
model.t2scfact=cmod.t2scfact; 

model.qlim=cmod.qlim;
model.qdof=cmod.qdof; 
model.qscfact=cmod.qscfact; 
model.dlim=cmod.dlim;

acc=zeros(ns,1); 

switch opt.cmcrit
    
    case 'sim'
        t2red=model.t2/cmod.t2lim; 
        qred=model.q/cmod.qlim; 
        dred=max(t2red, qred);
 
        
    case 'alt'
        
        t2red=model.t2/cmod.t2lim; 
        qred=model.q/cmod.qlim; 
        dred=sqrt((t2red.^2)+(qred.^2));

        
    case 'ci'
        
        t2red=model.t2/cmod.t2lim; 
        qred=model.q/cmod.qlim; 
        dred=t2red+qred;
        
    case 'dd'
        
        t2red=model.t2dof*model.t2/model.t2scfact; 
        qred=model.qdof*model.q/model.qscfact; 
        dred=t2red+qred;
        
end

        model.t2red=t2red; 
        model.qred=qred; 
        model.dred=dred; 
        
        acc(dred<=model.dlim)=1; 
        model.accepted=acc; 
        
        
        
    
