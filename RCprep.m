function [Xp, model]= RCprep(X, pstrat,prstr)
%RCprep Function for preprocessing data
%  The function RCprep sequentially preprocess data, according to the order
%  of pretreatments reported in the preprocessing structure prstr.
%  Preprocessing parameters can be based on the actual data, on a previous
%  training dataset %
%
%
%
% Author: Federico Marini
% Version: 21/04/2020



switch pstrat
    case 'model'
        [Xp,model]=Rcprepmod(X,prstr);
        
    case 'apply'
        
        [Xp,model]=Rcprepapply(X,prstr);
        
    case 'undo'
        
        Xp=Rcprepundo(X,prstr);
        model=prstr;
end


function [Xp,model]=Rcprepmod(X,prstr)
nprep=length(prstr);
model=prstr;
Xp=X;
[ns,nv]=size(Xp); 

for i=1:nprep
    switch prstr(i).type
        case 'none'
            Xp=Xp; 
            model(i).parameters=[];
            
        case 'mean'
            mx=mean(Xp); 
            model(i).parameters=mx; 
            Xp=Xp-repmat(mx,ns,1);
            
        case 'auto'
            mx=mean(Xp); 
            sx=std(Xp); 
            if isempty(prstr(i).settings) || ~strcmp(prstr(i).settings{1}, 'eps')
                Xp=Xp-repmat(mx,ns,1);
                Xp=Xp(:,sx~=0)./repmat(sx(sx~=0),ns,1);
            elseif ischar(prstr(i).settings{1}) && strcmp(prstr(i).settings{1}, 'eps')
                sx(sx==0)=eps; 
                Xp=(Xp-repmat(mx,ns,1))./repmat(sx,ns,1);
            elseif ~ischar(prstr(i).settings{1})
                sx(sx==0)=prstr(i).settings{1};
                Xp=(Xp-repmat(mx,ns,1))./repmat(sx,ns,1);
            end
            
            model(i).parameters={mx sx};           
            
        case 'pareto'
            
            mx=mean(Xp); 
            sx=sqrt(std(Xp)); 
            
            if isempty(prstr(i).settings) || ~strcmp(prstr(i).settings{1}, 'eps')
                Xp=Xp-repmat(mx,ns,1);
                Xp=Xp(:,sx~=0)./repmat(sx(sx~=0),ns,1);
            elseif ischar(prstr(i).settings{1}) && strcmp(prstr(i).settings{1}, 'eps')
                sx(sx==0)=eps; 
                Xp=(Xp-repmat(mx,ns,1))./repmat(sx,ns,1);
            elseif ~ischar(prstr(i).settings{1})
                sx(sx==0)=prstr(i).settings{1};
                Xp=(Xp-repmat(mx,ns,1))./repmat(sx,ns,1);
            end
            
            model(i).parameters={mx sx}; 
            
        case 'scalenc'
            
            sx=std(Xp); 
            if isempty(prstr(i).settings) || ~strcmp(prstr(i).settings{1}, 'eps')
                Xp=Xp(:,sx~=0)./repmat(sx(sx~=0),ns,1);
                Xp(:,sx==0)=0;
            elseif ischar(prstr(i).settings{1}) && strcmp(prstr(i).settings{1}, 'eps')
                sx(sx==0)=eps; 
                Xp=Xp./repmat(sx,ns,1);
            elseif ~ischar(prstr(i).settings{1})
                sx(sx==0)=prstr(i).settings{1};
                Xp=Xp./repmat(sx,ns,1);
            end
            
            model(i).parameters=sx;
            
            case 'paretonc'
            
            sx=sqrt(std(Xp)); 
            if isempty(prstr(i).settings) || ~strcmp(prstr(i).settings{1}, 'eps')
                Xp=Xp(:,sx~=0)./repmat(sx(sx~=0),ns,1);
                Xp(:,sx==0)=0;
            elseif ischar(prstr(i).settings{1}) && strcmp(prstr(i).settings{1}, 'eps')
                sx(sx==0)=eps; 
                Xp=Xp./repmat(sx,ns,1);
            elseif ~ischar(prstr(i).settings{1})
                sx(sx==0)=prstr(i).settings{1};
                Xp=Xp./repmat(sx,ns,1);
            end
            
            model(i).parameters=sx;
            
        case 'snv'
            mx=mean(Xp,2); 
            sx=std(Xp,[],2); 
            if isempty(prstr(i).settings) || ~strcmp(prstr(i).settings{1}, 'eps')
                Xp=Xp-repmat(mx,1,nv);
                Xp=Xp(sx~=0,:)./repmat(sx(sx~=0),1,nv);
            elseif ischar(prstr(i).settings{1}) && strcmp(prstr(i).settings{1}, 'eps')
                sx(sx==0)=eps; 
                Xp=(Xp-repmat(mx,1,nv))./repmat(sx,1,nv);
            elseif ~ischar(prstr(i).settings{1})
                sx(sx==0)=prstr(i).settings{1};
                Xp=(Xp-repmat(mx,1,nv))./repmat(sx,1,nv);
            end
            
            model(i).parameters={mx sx};
            
        case 'msc'
            if isempty(prstr(i).settings) || ~strcmp(prstr(i).settings{1}, 'mean')
                xref=mean(Xp); 
            elseif ischar(prstr(i).settings{1}) || ~strcmp(prstr(i).settings{1}, 'median')
                xref=median(Xp);
            elseif isnumeric(prstr(i).settings{1})
                xref=prstr(i).settings;
            end
            
            p=[ones(size(xref)); xref];
            C=Xp*pinv(p);
            Xp=(Xp-repmat(C(:,1), 1, 3112))./repmat(C(:,2), 1, 3112);
            
            model(i).parameters={xref, C(:,1), C(:,2)};
            
        case 'detrend'
            if isempty(prstr(i).settings)
                model(i).settings{1}=0; 
            end
            
            p=zeros(model(i).settings{1}+1,nv);
            p(1,:)=ones(1,nv);
            if model(i).settings{1}>0
                if length(model(i).settings)>1 && ~isempty(model(i).settings{2})
                    xsc=model(i).settings{2};
                else
                    xsc=1:nv; 
                end
                for j=1:model(i).settings{1}
                    p(j+1,:)=xsc.^j;
                end
            end
            C=Xp*pinv(p);
            Xp=Xp-C*p; 
            model(i).parameters={p,C};
            
        case 'norm'
            if isempty(prstr(i).settings)
                model(i).settings{1}='L2'; 
            end
            
            switch model(i).settings{1}
                case 'L1'
                    C=sum(Xp,2); 
                case 'L2'
                    C=sqrt(sum(Xp.^2,2)); 
                case 'max'
                    C=max(Xp,[],2); 
            end
            
            Xp=Xp./repmat(C,1,nv); 
            model(i).parameters=C;
            
        case 'pqn'
            if isempty(prstr(i).settings)
                model(i).settings={'median'}; 
                xref=[];
            elseif ~isempty(prstr(i).settings) && isnumeric(prstr(i).settings{1})
                xref=prstr(i).settings{1};
            else
                xref=[]; 
                if length(prstr(i).settings)>1 && isnumeric(prstr(i).settings{2})
                    sels=prstr(i).settings{2}; 
                else
                    sels=1:ns;
                end
            end
            
            if ~isempty(xref)
                switch model(i).settings
                    case 'median'
                        xref=median(Xp(sels,:)); 
                    case 'mean'
                        xref=mean(Xp(sels,:)); 
                    case 'max'
                        xref=max(Xp(sels,:));
                        case 'min'
                        xref=min(Xp(sels,:));
                end
            end
            %calculation of coefficients
            C=median(Xp./repmat(xref, ns,1),2);
            Xp=Xp./repmat(C, 1, nv);
            
            model(i).parameters={xref,C};

            
        case 'vsn'
            if isempty(prstr(i).settings)
                [Xp, vsnres]=vsn(Xp);
            else
                [Xp, vsnres]=vsn(Xp,prstr(i).settings{1});
            end
            model(i).parameters=vsnres; 
            
        
        case 'savgol'
            if isempty(prstr(i).settings)
                model(i).settings={7,2,0};    
            end
            Xp=RCsavgol(Xp, model(i).settings{1}, model(i).settings{2}, model(i).settings{3});
            
        case 'AsLS'
            if isempty(prstr(i).settings)
                model(i).settings={10^5,0.01};
            end
            
            bline=zeros(size(Xp));
            for j=1:ns
                [Xp(j,:), bline(j,:)]=baseline_als(Xp(j,:),model(i).settings{1} ,model(i).settings{1}, 10);
            end
            model(i).parameters=bline;
            
        case 'blocksc'
            sb=norm(Xp, 'fro');
            Xp=Xp/sb; 
            model(i).parameters=sb;
            
        case 'TtoA'
            if isempty(prstr(i).settings) || ~strcmp(prstr(i).settings{1}, '%')
                Xp=-log10(Xp);
            else
                Xp=-log10(0.01*Xp);
                
            end
            
            
            
    end
    
end


function [Xp,model]=Rcprepapply(X,prstr)
nprep=length(prstr);
model=prstr;
Xp=X;
[ns,nv]=size(Xp); 

for i=1:nprep
    switch prstr(i).type
        case 'none'
            Xp=Xp; 
            
        case 'mean'
            mx=prstr(i).parameters; 
            Xp=Xp-repmat(mx,ns,1);
            
        case 'auto'
            mx=model(i).parameters{1}; 
            sx=model(i).parameters{2}; 
            if isempty(prstr(i).settings)
                Xp=Xp-repmat(mx,ns,1);
                Xp=Xp(:,sx~=0)./repmat(sx(sx~=0),ns,1);
            else
                Xp=(Xp-repmat(mx,ns,1))./repmat(sx,ns,1);
            end
                      
            
        case 'pareto'
            
            mx=model(i).parameters{1}; 
            sx=model(i).parameters{2}; 
            if isempty(prstr(i).settings)
                Xp=Xp-repmat(mx,ns,1);
                Xp=Xp(:,sx~=0)./repmat(sx(sx~=0),ns,1);
            else
                Xp=(Xp-repmat(mx,ns,1))./repmat(sx,ns,1);
            end
            
        case 'scalenc'
            
            sx=model(i).parameters; 
            if isempty(prstr(i).settings)
                Xp=Xp(:,sx~=0)./repmat(sx(sx~=0),ns,1);
                Xp(:,sx==0)=0;
            else
                Xp=Xp./repmat(sx,ns,1);
            end
            
            
            
            case 'paretonc'
            
            sx=model(i).parameters; 
            if isempty(prstr(i).settings)
                Xp=Xp(:,sx~=0)./repmat(sx(sx~=0),ns,1);
                Xp(:,sx==0)=0;
            else
                Xp=Xp./repmat(sx,ns,1);
            end
            
            
        case 'snv'
            mx=mean(Xp,2); 
            sx=std(Xp,[],2); 
            if isempty(prstr(i).settings) || ~strcmp(prstr(i).settings{1}, 'eps')
                Xp=Xp-repmat(mx,1,nv);
                Xp=Xp(sx~=0,:)./repmat(sx(sx~=0),1,nv);
            elseif ischar(prstr(i).settings{1}) && strcmp(prstr(i).settings{1}, 'eps')
                sx(sx==0)=eps; 
                Xp=(Xp-repmat(mx,1,nv))./repmat(sx,1,nv);
            elseif ~ischar(prstr(i).settings{1})
                sx(sx==0)=prstr(i).settings{1};
                Xp=(Xp-repmat(mx,1,nv))./repmat(sx,1,nv);
            end
            
            model(i).parameters={mx sx};
            
        case 'msc'
             
            p=[ones(size(xref)); prstr(i).parameters{1}];
            C=Xp*pinv(p);
            Xp=(Xp-repmat(C(:,1), 1, 3112))./repmat(C(:,2), 1, 3112);
            
            model(i).parameters(2:3)={C(:,1), C(:,2)};
            
        case 'detrend'
            C=Xp*pinv(prstr(i).parameters{1});
            Xp=Xp-C*prstr(i).parameters{1}; 
            model(i).parameters{2}=C;
            
        case 'norm'
            switch prstr(i).settings{1}
                case 'L1'
                    C=sum(Xp,2); 
                case 'L2'
                    C=sqrt(sum(Xp.^2,2)); 
                case 'max'
                    C=max(Xp,[],2); 
            end
            
            Xp=Xp./repmat(C,1,nv); 
            model(i).parameters=C;
            
        case 'pqn'
            
            C=median(Xp./repmat(prstr(i).parameters{1}, ns,1),2);
            Xp=Xp./repmat(C, 1, nv);
            
            model(i).parameters{2}=C;

        case 'vsn'
            Xp=vsn(Xp, prstr(i).parameters); 
            
        case 'savgol'
            Xp=RCsavgol(Xp, prstr(i).settings{1}, prstr(i).settings{2}, prstr(i).settings{3}); 
            
            
            
        case 'AsLS'
            
            bline=zeros(size(Xp));
            for j=1:ns
                [Xp(j,:), bline(j,:)]=baseline_als(Xp(j,:),prstr(i).settings{1} ,prstr(i).settings{1}, 10);
            end
            model(i).parameters=bline;
            
        case 'blocksc'
            Xp=Xp/prstr(i).parameters; 
            
        case 'TtoA'
            if isempty(prstr(i).settings) || ~strcmp(prstr(i).settings{1}, '%')
                Xp=-log10(Xp);
            else
                Xp=-log10(0.01*Xp);
                
            end
            
            
            
    end
    
end

function Xp=Rcprepundo(X,prstr)
nprep=length(prstr);
Xp=X;
[ns,nv]=size(Xp); 

for i=nprep:-1:1
    switch prstr(i).type
        case 'none'
            Xp=Xp; 
            
        case 'mean'
            mx=prstr(i).parameters; 
            Xp=Xp+repmat(mx,ns,1);
            
        case 'auto'
            mx=prstr(i).parameters{1}; 
            sx=prstr(i).parameters{2}; 
            Xp=(Xp.*repmat(sx,ns,1))+repmat(mx,ns,1);
            
        case 'pareto'
            
            mx=prstr(i).parameters{1}; 
            sx=prstr(i).parameters{2}; 
            Xp=(Xp.*repmat(sx,ns,1))+repmat(mx,ns,1);
            
        case 'scalenc'
            
            sx=prstr(i).parameters; 
            Xp=Xp.*repmat(sx,ns,1);
            
            
            
            case 'paretonc'
            
            sx=prstr(i).parameters; 
            Xp=Xp.*repmat(sx,ns,1);
            
            
        case 'snv'
            mx=prstr(i).parameters{1}; 
            sx=prstr(i).parameters{2}; 
            
            Xp=(Xp.*repmat(sx,1,nv))+repmat(mx,1,nv);
            
            
        case 'msc'
             
            p=[ones(size(xref)); prstr(i).parameters{1}];
            C=Xp*pinv(p);
            Xp=(Xp.*repmat(prstr(i).parameters{3}, 1, 3112))+repmat(prstr(i).parameters{2}, 1, 3112);
            
            
        case 'detrend'
            Xp=Xp+prstr(i).parameters{2}*prstr(i).parameters{1}; 
            
            
        case 'norm'
            Xp=Xp.*repmat(prstr(i).parameters,1,nv); 
            
        case 'pqn'
            
           Xp=Xp.*repmat(prstr(i).parameters{2}, 1, nv);
            
            
        case 'savgol'
            Xp=Xp; 
            
            
            
        case 'AsLS'
            
            Xp=Xp+prstr(i).parameters;
            
        case 'blocksc'
            Xp=Xp*prstr(i).parameters; 
            
        case 'TtoA'
            if isempty(prstr(i).settings) || ~strcmp(prstr(i).settings{1}, '%')
                Xp=10.^(-Xp);
            else
                Xp=10.^(-100*Xp);
                
            end
            
            
            
    end
    
end

function ynew=RCsavgol(varargin)
%RCsavgol Function for Savitzky-Golay smoothing and differentiation
%
% Author: Federico Marini
% 
if nargin==4
    y=varargin{1}; 
    width=varargin{2}; 
    order=varargin{3}; 
    deriv=varargin{4};
elseif nargin==3
    y=varargin{1}; 
    width=varargin{2}; 
    order=varargin{3}; 
    deriv=0;
elseif nargin==2
    y=varargin{1}; 
    width=varargin{2}; 
    order=2; 
    deriv=0;    
elseif nargin==1
    y=varargin{1}; 
    width=7; 
    order=2; 
    deriv=0;
end

[r,c] = size(y);
ynew=y;

polord=min([max(0,round(order)),width-1]);
der = min(max(0,round(deriv)),polord);


hpoints = (width-1)/2;

% Building the independent block and its pseudoinverse
X = repmat((-hpoints:hpoints)',1,1+polord).^repmat((0:polord), width,1);
A = ((X'*X)\X')';

% Calcuating the SG filter/derivative for almost all data
for i=hpoints+2:c-hpoints-1                              
  ynew(:,i) = y(:,i-hpoints:i+hpoints)*(prod(1:der)*A(:,der+1));      
end

% SG smoothing/derivative for the tails
A = [y(:,1:width); y(:,c-width+1:c)]*A;  
for i=1:der
  A = A(:,2:polord+2-i)*diag(1:polord+1-i); % or its d'th derivative
end
ynew(:,1:hpoints+1) = A(1:r,:)*X(1:hpoints+1,1:1+polord-der)'; 
ynew(:,c-hpoints:c) = A(r+(1:r),:)*X(hpoints+1:width,1:1+polord-der)';


function [ycorr,z] = baseline_als(y, lambda, p, maxit)
% Estimate baseline with asymmetric least squares
% I/O: [ycorr,z] = baseline_als(y, lambda, p, maxit)
% common values for lambda are 10^5-10^8
% while for p 0.005-0.01
y=y';

m = length(y);
D = diff(speye(m), 2);
w = ones(m, 1);
for it = 1:maxit
W = spdiags(w, 0, m, m);
C = chol(W + lambda * (D' * D));
z = C \ (C' \ (w .* y));
w = p * (y > z) + (1 - p) * (y < z);
end
z=z';
ycorr=y'-z;

        