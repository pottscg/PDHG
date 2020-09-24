%% Morup's S (A) update function
% this method uses projected gradient vs PDHG

function [S,SSE,muS,SSt]=Supdate(S,XCtX,CtXtXC,muS,SST,SSE,niter)
    
    [noc,J]=size(S);
    e=ones(noc,1);
    for k=1:niter
        SSE_old=SSE;
        g=(CtXtXC*S-XCtX)/(SST/J);         
        g=g-e*sum(g.*S,1);
        stop=0;
        Sold=S;
        while ~stop
            S=Sold-g*muS;
            S(S<0)=0;
            S=S./(e*sum(S,1));            
            SSt=S*S';          
            SSE=SST-2*sum(sum(XCtX.*S))+sum(sum(CtXtXC.*SSt));            
            if SSE<=SSE_old*(1+1e-9)
                muS=muS*1.2;
                stop=1;
            else
                muS=muS/2;
            end
        end
    end