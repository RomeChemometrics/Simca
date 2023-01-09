function model=RCsimcapcapred(X, simpcmod)

% Preprocessing routine
Xc=RCprep(X, 'apply', simpcmod.pretX.preppars); 
[ns,nv]=size(Xc); 

T=Xc*simpcmod.loadings; 
model=simpcmod; 

model.scores=T; 

% Calculation of Q contribution and of Q statistics
if ns>nv
    qcon=Xc*(eye(nv)-simpcmod.loadings*simpcmod.loadings');
    q=sum(qcon.^2,2);
    qcon=sign(qcon).*qcon.^2;
else
    qcon=Xc-T*simpcmod.loadings'; 
    q=sum(qcon.^2,2);
    qcon=sign(qcon).*qcon.^2;
    
end

% Calculation of T2 and its contribution
t2con = T/(diag(sqrt(simpcmod.eigs)))*simpcmod.loadings';
t2=sum(t2con.^2,2);
t2con=sign(t2con).*t2con.^2;

model.t2=t2; 
model.t2con=t2con; 
model.q=q; 
model.qcon=qcon; 
