function [C,SSE,muC,mualpha,CtXtXC,XC]=Cupdate(X,XSt,XC,SSt,C,delta,muC,mualpha,SST,SSE,niter)
                                       
    [J,noc]=size(C);
    if nargin<12
        niter=1;
    end   
    if delta~=0
        alphaC=sum(C);    
        C=C*diag(1./alphaC);
    end
    e=ones(J,1);
    XtXSt=X'*XSt;
    for k=1:niter
        
        % Update C        
        SSE_old=SSE;        
        g=(X'*(XC*SSt)-XtXSt)/SST;       
        if delta~=0
            g=g*diag(alphaC);
        end
        g=g-e*sum(g.*C,1);        
        stop=0;
        Cold=C;
        while ~stop
            C=Cold-muC*g;
            C(C<0)=0;
            
            nC=sum(C,1)+eps;            
            C=C*sparse(1:noc,1:noc,1./nC);
            if delta~=0
                Ct=sparse(C*diag(alphaC));
            else
                Ct=sparse(C);
            end
            XC=X*Ct;
            CtXtXC=XC'*XC;
            SSE=SST-2*sum(sum(XC.*XSt))+sum(sum(CtXtXC.*SSt));                                    
            if SSE<=SSE_old*(1+1e-9)
                muC=muC*1.2;
                stop=1;
            else
                muC=muC/2;
            end
        end        
        
        % Update alphaC        
        SSE_old=SSE;
        if delta~=0                                                           
            g=(diag(CtXtXC*SSt)'./alphaC-sum(C.*XtXSt,1))/(SST*J);                       
            stop=0;
            alphaCold=alphaC;
            while ~stop
                alphaC=alphaCold-mualpha*g;
                alphaC(alphaC<1-delta)=1-delta;
                alphaC(alphaC>1+delta)=1+delta;                            
                XCt=XC*diag(alphaC./alphaCold);
                CtXtXC=XCt'*XCt;            
                SSE=SST-2*sum(sum(XCt.*XSt))+sum(sum(CtXtXC.*SSt));                
                if SSE<=SSE_old*(1+1e-9)  
                    mualpha=mualpha*1.2;
                    stop=1;
                    XC=XCt;
                else
                    mualpha=mualpha/2;
                end
            end  
        end
    end
    if delta~=0
        C=C*diag(alphaC);
    end

