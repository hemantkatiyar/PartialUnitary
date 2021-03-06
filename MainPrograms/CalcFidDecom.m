% Author  : Dr. Hemant Katiyar
% Email   : hkatiyar@uwaterloo.ca
% Website : 

% Description :
%

function [Fid,Uf] = CalcFidDecom(x,nSpin,nSec,VarPerSec,DelayControl,Iz,Had,Hint,Utarg,D)

% Function Definition
CalcFid = @(x,y,z) abs(trace(x*y'))/z;
sTarg=size(Utarg);
D1=sTarg(2);

Uf  = cell(1,2*nSec+2);
if nSec~=0
    for n = 1:nSec
        Ievo = n*VarPerSec;
        Ilast = (n-1)*VarPerSec;
        
        PhIz = sum(repmat(x(Ilast+1:Ilast+nSpin),D,1).*Iz,2);
        AIz  = sum(repmat(x(Ilast+nSpin+1:Ilast+2*nSpin),D,1).*Iz,2);
    
        Rxy  = (exp(-1i*PhIz).*exp(-1i*PhIz)').*(Had*(exp(-1i*AIz).*Had'));
        Revo = exp(-1i*abs(x(Ievo))*DelayControl/pi*Hint);
    
        if n==1
            Uf{1,2*n-1} = Rxy;
        else
            Uf{1,2*n-1} = Rxy*Uf{1,2*n-2};
        end
        Uf{1,2*n}   = Revo.*Uf{1,2*n-1};    
    end
    
    Ilast=VarPerSec*nSec;
    PhIz = sum(repmat(x(Ilast+1:Ilast+nSpin),D,1).*Iz,2);
    AIz  = sum(repmat(x(Ilast+nSpin+1:Ilast+2*nSpin),D,1).*Iz,2);
    AzIz = sum(repmat(x(Ilast+2*nSpin+1:Ilast+3*nSpin),D,1).*Iz,2);
    
    Rxy  = (exp(-1i*PhIz).*exp(-1i*PhIz)').*(Had*(exp(-1i*AIz).*Had'));
    Rz   = exp(-1i*AzIz);
    
    Uf{1,2*nSec+1}=Rxy*Uf{1,2*nSec};
    Uf{1,2*nSec+2}=Rz.*Uf{1,2*nSec+1};
    
    Fid = CalcFid(Uf{1,end}(1:sTarg(1),1:sTarg(2)),Utarg,D1);
else % Only 2 Unitaries Needed
    Ilast=0;
    PhIz = sum(repmat(x(Ilast+1:Ilast+nSpin),D,1).*Iz,2);
    AIz  = sum(repmat(x(Ilast+nSpin+1:Ilast+2*nSpin),D,1).*Iz,2);
    AzIz = sum(repmat(x(Ilast+2*nSpin+1:Ilast+3*nSpin),D,1).*Iz,2);
    
    Rxy  = (exp(-1i*PhIz).*exp(-1i*PhIz)').*(Had*(exp(-1i*AIz).*Had'));
    Rz   = exp(-1i*AzIz);
    
    Uf{1,1}=Rxy;
    Uf{1,2}=Rz.*Uf{1,1};
    Fid = CalcFid(Uf{1,end}(1:sTarg(1),1:sTarg(2)),Utarg,D1);
end


