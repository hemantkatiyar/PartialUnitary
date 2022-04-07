% Author  : Dr. Hemant Katiyar
% Email   : hkatiyar@uwaterloo.ca
% Website : 

% Description :
% 
function Grad = CalcGradDecom(x,nSpin,nSec,VarPerSec,DelayControl,Utarg,Ix,Iy,Iz,Hint,Uf,D)

Ub=cell(1,2*nSec+2);
for kk=1:2*nSec+2
    Ub{1,kk}=Uf{1,end}*Uf{1,kk}';
end

EvalGrad = @(x,y,z,a,p) (2*real(trace(x(1:p(1),1:p(2))*y')*z)/a^2);

p=size(Utarg);
D=p(2);

trxjpj = trace(Uf{1,end}(1:p(1),1:p(2))'*Utarg);
Grad= zeros(size(x));
if nSec~=0
    for n=1:nSec
        Ievo = n*VarPerSec;
        Ilast = (n-1)*VarPerSec;
        
        for j = 1:nSpin
            if n==1
                Grad(Ilast+j) = ...
                EvalGrad(-1i*(Ub{1,2*n-1}*(Iz(:,j).*Uf{1,2*n-1})-Uf{1,end}.*transpose(Iz(:,j))),Utarg,trxjpj,D,p);
            else
                Grad(Ilast+j) = ...
                EvalGrad(-1i*(Ub{1,2*n-1}*(Iz(:,j).*Uf{1,2*n-1})-Ub{1,2*n-2}*(Iz(:,j).*Uf{1,2*n-2})),Utarg,trxjpj,D,p);
            end
            Grad(Ilast+nSpin+j) = ...
            EvalGrad(-1i*Ub{1,2*n-1}*(cos(x(Ilast+j))*Ix(:,:,j)+sin(x(Ilast+j))*Iy(:,:,j))*Uf{1,2*n-1},Utarg,trxjpj,D,p);
        end
        Grad(Ievo) = ...
        EvalGrad(-1i*Ub{1,2*n}*(((x(Ievo)/abs(x(Ievo)))*DelayControl/pi*Hint).*Uf{1,2*n}),Utarg,trxjpj,D,p);
    end
    
    Ilast=VarPerSec*nSec;
    for j = 1:nSpin
        Grad(Ilast+j) = ...
        EvalGrad(-1i*(Ub{1,end-1}*(Iz(:,j).*Uf{1,end-1})-Ub{1,end-2}*(Iz(:,j).*Uf{1,end-2})),Utarg,trxjpj,D,p);
        Grad(Ilast+nSpin+j) = ...
        EvalGrad(-1i*Ub{1,end-1}*(cos(x(Ilast+j))*Ix(:,:,j)+sin(x(Ilast+j))*Iy(:,:,j))*Uf{1,end-1},Utarg,trxjpj,D,p); 
        Grad(Ilast+2*nSpin+j) = EvalGrad(-1i*(Iz(:,j).*Uf{1,end}),Utarg,trxjpj,D,p);
    end
else
    Ilast=0;
    for j = 1:nSpin
        Grad(Ilast+j) = ...
        EvalGrad(-1i*(Ub{1,end-1}*(Iz(:,j).*Uf{1,end-1})-Uf{1,end}.*transpose(Iz(:,j))),Utarg,trxjpj,D,p);
        Grad(Ilast+nSpin+j) = ...
        EvalGrad(-1i*Ub{1,end-1}*(cos(x(Ilast+j))*Ix(:,:,j)+sin(x(Ilast+j))*Iy(:,:,j))*Uf{1,end-1},Utarg,trxjpj,D,p); 
        Grad(Ilast+2*nSpin+j) = EvalGrad(-1i*(Iz(:,j).*Uf{1,end}),Utarg,trxjpj,D,p);
    end
end


% Ub{1} =  Uf{2}*Uf{1}' = U{2}
% Ub{2} =  Uf{2}*Uf{2}' = Id
% Uf{1} = U1
% Uf{2} = U2*U1


% Ub1 Iz Uf1 - Ub0 Iz Uf0
% U2  Iz U1  - U2U1  Iz Id
