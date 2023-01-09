function model=RCsimcapcamod(X, nF, pX)

% Preprocessing routine

nP=length(pX);
pstr=struct('type', [], 'settings', []);
for j=1:nP
    
    if ischar(pX{j})
        pstr(j).type=pX{j};
    elseif iscell(pX{j})
        pstr(j).type=pX{j}{1};
        pstr(j).settings=pX{j}(2:end);
    end
end

[Xc, ppars]=RCprep(X, 'model', pstr);

[ns,nv]=size(Xc); 

nF=min(nF, rank(Xc)); 

%Computes nF factors PCA model
if ns>nv
    [~,S,P]=svd(Xc'*Xc, 'econ');
    P=P(:,1:nF);
    T=Xc*P;    
    S=diag(S)./(ns-1);
else
    [u,S,P]=svd(Xc,'econ');
    P=P(:,1:nF);
    T=u(:,1:nF)*S(1:nF,1:nF);
    S=diag(S.^2)./(ns-1);
    
end

model.scores=T; 
model.loadings=P; 

lam=S; % All eigenvalues; 
model.eigs=S(1:nF); %Only significant eigenvalues
model.nPC=nF; 
%variance explained by each PC
ev=100*lam/sum(lam);
%cumulative variance explained by the first k Pcs
cv=cumsum(ev);

model.expl_var=ev(1:nF);
model.cum_var=cv(1:nF);


% Calculation of Q contribution and of Q statistics
if ns>nv
    qcon=Xc*(eye(nv)-P*P');
    q=sum(qcon.^2,2);
    qcon=sign(qcon).*qcon.^2;
else
    qcon=Xc-T*P'; 
    q=sum(qcon.^2,2);
    qcon=sign(qcon).*qcon.^2;
    
end

% Calculation of T2 and its contribution
t2con = T/(diag(sqrt(model.eigs)))*P';
t2=sum(t2con.^2,2);
t2con=sign(t2con).*t2con.^2;

model.t2=t2; 
model.t2con=t2con; 
model.q=q; 
model.qcon=qcon; 

model.pretX.ptype=pX; 
model.pretX.details=pstr; 
model.pretX.preppars=ppars; 
model.detail.eigs=lam; 
model.detail.expl_var=ev; 
model.detail.cum_var=cv; 


