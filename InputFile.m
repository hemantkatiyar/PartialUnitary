;Par=[];
% Name of molecule file located in MoleculeFiles folder
Par.MoleculeFile='MolCrotonic7';
% Par.MoleculeFile='MolCrotonic7a';
% Par.MoleculeFile='MolCrotonic';
% Par.MoleculeFile='SvQubitError
% Par.MoleculeFile='SvQubit';
% Par.MoleculeFile='TwQubit';
% cnot = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0];
% Par.Utarg = kron(cnot,eye(2^2));


% Number of sections in decomposition
% Par.nSec = 10;

% Max value of sum of all delay times for the guess (in ms)
Par.MaxIniDelay = 100;

% Parmeter to control t
Par.DelayControl = 0.001;
% Initial step size(usually left alone)
Par.init_stepsize = 1;

% If GuessFlag=0, then nohe delay times(usually left alone)
% Par. guess is used, if its 1 load the GuessFileName

% Par.GuessFlag=0;
% Par.GuessFileName='CR7_PPS_Encoding10';
% If algorithm stops prematurely, should we keep trying? if yes, then 1
Par.KeepTrying=0;

% Save the results wiht the following name in SaveOutputs\SaveOutputsDecom

% Define Unitary that needs decomposing


% [Ix,Iy,Iz,~,~,~,~] = prodopSparse([1/2],[4]);
% Par.Utarg = expm(1i*pi/2*Iy(:,:,1))*(exp(-1i*pi*Iz(:,1).*Iz(:,2)).*expm(-1i*pi/2*Ix(:,:,1))); Par.SaveFileName='UZ1Z2s';
% Par.Utarg = expm(1i*pi/2*Iy(:,:,2))*(exp(-1i*pi*Iz(:,2).*Iz(:,3)).*expm(-1i*pi/2*Ix(:,:,2))); 
% Par.Utarg = expm(1i*pi/4*Iy(:,:,3))*(exp(-1i*pi*Iz(:,3).*Iz(:,4)).*expm(-1i*pi/4*Ix(:,:,3)));
% Par.Utarg = expm(1i*pi/4*Iy(:,:,2))*(exp(-1i*pi*Iz(:,2).*Iz(:,3)).*expm(-1i*pi/4*Ix(:,:,2))); Par.SaveFileName='UZ2Z3s';
% Par.Utarg = expm(1i*pi/4*Iy(:,:,1))*(exp(-1i*pi*Iz(:,1).*Iz(:,2)).*expm(-1i*pi/4*Ix(:,:,1))); Par.SaveFileName='UZ1Z2_45s';

% Par.SaveFileName='U1EvoM';
% [Ix,Iy,Iz,~,~,~,~] = prodopSparse([1/2],[4]);
% Par.Utarg = expm(1i*pi/2*Iy(:,:,1))*expm(-1i*pi*(Iy(:,:,1)+Iy(:,:,2)))*(exp(-1i*pi*Iz(:,1).*Iz(:,2)).*expm(-1i*pi/2*Ix(:,:,1)));
% Par.Utarg = expm(1i*pi/2*Iy(:,:,2))*expm(-1i*pi*(Iy(:,:,2)+Iy(:,:,3)))*(exp(-1i*pi*Iz(:,2).*Iz(:,3)).*expm(-1i*pi/2*Ix(:,:,2)));
% Par.Utarg = expm(-1i*pi/4*Iy(:,:,3))*expm(-1i*pi*(Iy(:,:,3)+Iy(:,:,4)))*(exp(-1i*pi*Iz(:,3).*Iz(:,4)).*expm(-1i*pi/4*Ix(:,:,3)));
% Par.Utarg = expm(-1i*pi/4*Iy(:,:,2))*expm(-1i*pi*(Iy(:,:,2)+Iy(:,:,3)))*(exp(-1i*pi*Iz(:,2).*Iz(:,3)).*expm(-1i*pi/4*Ix(:,:,2)));
% Par.Utarg = expm(-1i*pi/4*Iy(:,:,1))*expm(-1i*pi*(Iy(:,:,1)+Iy(:,:,2)))*(exp(-1i*pi*Iz(:,1).*Iz(:,2)).*expm(-1i*pi/4*Ix(:,:,1)));

% [Ix,Iy,Iz,~,~,~,~] = prodopSparse([1/2],[7]);
% % Par.Utarg = ((exp(-1i*pi*Iz(:,1).*Iz(:,2)).*exp(-1i*pi*Iz(:,3).*Iz(:,2))).*expm(-1i*pi/2*Iy(:,:,2))); Par.SaveFileName='UZ1Z2s';
% 
% Par.Utarg=expm(-1i*pi/2*Iy(:,:,2))*(exp(-1i*pi*Iz(:,2).*Iz(:,3)).*(exp(-1i*pi*Iz(:,2).*Iz(:,1)).*(expm(-1i*pi/2*Iy(:,:,2))*(expm(1i*pi/2*Iy(:,:,7))*(exp(-1i*pi*Iz(:,2).*Iz(:,7)).*(exp(-1i*pi*Iz(:,4).*Iz(:,7)).*(exp(-1i*pi*Iz(:,5).*Iz(:,7)).*(exp(-1i*pi*Iz(:,6).*Iz(:,7)).*expm(-1i*pi/2*Iy(:,:,7))))))))));
% A1=Par.Utarg*diag(Iz(:,7))*Par.Utarg';
% A2=diag(Iz(:,1))*diag(Iz(:,2))*diag(Iz(:,3).*Iz(:,4).*Iz(:,5).*Iz(:,6).*Iz(:,7));
% trace(A1*A2)/sqrt(trace(A1^2)*trace(A2^2))
% Par.SaveFileName='UEncoding8test';
% Par.Utarg=expm(1i*pi/2*(Iy(:,:,2)+Iy(:,:,5)))*(exp(-1i*pi*Iz(:,5).*Iz(:,7)).*(exp(-1i*pi*Iz(:,2).*Iz(:,7)).*(expm(-1i*pi/2*(Ix(:,:,2)+Ix(:,:,5)))*expm(1i*pi/2*(Iy(:,:,1)+Ix(:,:,3)+Ix(:,:,4)+Ix(:,:,6)))*(exp(-1i*pi*Iz(:,5).*Iz(:,6)).*(exp(-1i*pi*Iz(:,4).*Iz(:,5)).*(exp(-1i*pi*Iz(:,3).*Iz(:,2)).*(exp(-1i*pi*Iz(:,1).*Iz(:,2)).*expm(-1i*pi/2*(Ix(:,:,1)+Ix(:,:,3)+Ix(:,:,4)+Ix(:,:,6))))))))));
% Par.Utarg = (expm(1i*pi/2*(Iy(:,:,5)))*(exp(-1i*pi*Iz(:,5).*Iz(:,7)).*(expm(1i*pi/2*(Iy(:,:,2)-Ix(:,:,5)))*(exp(-1i*pi*Iz(:,2).*Iz(:,7)).*(expm(1i*pi/2*(Iy(:,:,6)-Ix(:,:,2)))*(exp(-1i*pi*Iz(:,5).*Iz(:,6)).*(expm(1i*pi/2*(Iy(:,:,4)-Ix(:,:,6)))*(exp(-1i*pi*Iz(:,4).*Iz(:,5)).*(expm(1i*pi/2*(Iy(:,:,3)-Ix(:,:,4)))*(exp(-1i*pi*Iz(:,2).*Iz(:,3)).*(expm(1i*pi/2*(Iy(:,:,1)-Ix(:,:,3)))*(exp(-1i*pi*Iz(:,1).*Iz(:,2)).*expm(-1i*pi/2*Ix(:,:,1))))))))))))));
% Rin=zeros(2^7); Rin(128,1)=1; Rin(1,128)=1;
% A1=Par.Utarg*Rin*Par.Utarg';
% Rout=zeros(2^6);Rout(1,1)=1;
% A2=kron(Rout,[0 1;1 0]);
% trace(A1*A2)/sqrt(trace(A1^2)*trace(A2^2))
% sf
% Par.SaveFileName='UDecoding7';

% trace(A1*A2)/sqrt(trace(A1^2)*trace(A2^2))
% returndsf
% U=eye(2^7);
% U(128,128)=0; U(127,127)=0;U(128,127)=1; U(127,128)=1;
% x=[0 1; 1 0];
% Par.Utarg=kron(x,kron(x,kron(x,kron(x,kron(x,kron(x,x))))))*U;
% U=eye(2^3);
% U(8,8)=0; U(7,7)=0;U(8,7)=1; U(7,8)=1;
% Par.Utarg=kron(U,eye(2^4));
% % trace(kron(U,eye(2^4)))/2^7
% Par.SaveFileName='U7aSV';
%---------------------------------------------------------------------
% O=[1;0]; l=[0;1];
% A=O*O';  B=O*l'; C=l*O';  D=l*l';
% X=[0 1; 1 0];
% % 271 273 274 276
% % 127 237 247 267
% % U1 = kron(kron(kron(X,A),eye(2^4)),A)+kron(kron(kron(eye(2),A),eye(2^4)),D)+kron(kron(kron(eye(2),D),eye(2^4)),A)+kron(kron(kron(eye(2),D),eye(2^4)),D);
% U1 = kron(kron(kron(X,A),eye(2^4)),A) + kron(kron(kron(kron(eye(2),A),X),eye(2^3)),D)...
%     +kron(kron(kron(kron(kron(eye(2),D),eye(2)),X),eye(2^2)),A)+...
%      kron(kron(kron(kron(eye(2),D),eye(2^3)),X),D);
% % U1 = kron(kron(kron(X,A),eye(2^4)),A) + ...
%      % kron(kron(kron(eye(2),A),eye(2^4)),D) + ...
%      % kron(kron(kron(eye(2),D),eye(2^4)),A) + ...
%      % kron(kron(kron(eye(2),D),eye(2^4)),D);
% 
% % max(U1*U1'-eye(2^7))
% Par.Utarg=U1;
% Par.SaveFileName='UStep1SVc';
%---------------------------------------------------------------------
% O=[1;0]; l=[0;1];X=[0 1; 1 0]; I=eye(2); II=eye(2^2); III=eye(2^3);
% ll=kron(l,l)*kron(l',l');
% Toffoli236 = kron(kron(I,kron(II-ll,kron(II,I))),I)+kron(kron(I,kron(ll,kron(II,X))),I);
% Par.Utarg=Toffoli236;
% Par.SaveFileName='CR7_Toffoli236';
% Par.GuessFlag=1;
% Par.GuessFileName='CR7_Toffoli236_ss';
% Par.nSec = 11;
%---------------------------------------------------------------------
% O=[1;0]; l=[0;1];X=[0 1; 1 0]; I=eye(2); II=eye(2^2); III=eye(2^3);
% ll=kron(l,l)*kron(l',l'); OO = kron(O,O)*kron(O',O'); IIm = II-ll-OO;
% Toffoli236_231 = kron(kron(I,kron(IIm,kron(II,I))),I)+kron(kron(X,kron(OO,kron(II,I))),I)+kron(kron(I,kron(ll,kron(II,X))),I);
% Par.Utarg=Toffoli236_231;
% Par.SaveFileName='CR7_Toffoli236_231';
% Par.GuessFlag=1;
% Par.GuessFileName='CR7_Toffoli236_231_ss';
% Par.nSec = 15;
%---------------------------------------------------------------------
% O=[1;0]; l=[0;1];X=[0 1; 1 0]; I=eye(2); II=eye(2^2); III=eye(2^3);
% Ol=kron(O,l)*kron(O',l'); lO = kron(l,O)*kron(l',O'); IIm = II-Ol-lO;
% Toffoli235_234 = kron(kron(I,kron(IIm,kron(II,I))),I)+kron(kron(I,kron(Ol,kron(X,I))),II)+kron(kron(I,kron(lO,kron(I,X))),II);
% Par.Utarg=Toffoli235_234;
% Par.SaveFileName='CR7_Toffoli235_234_17sec';
% Par.GuessFlag=1;
% Par.GuessFileName='CR7_Toffoli235_23417';
% Par.nSec = 17;
%---------------------------------------------------------------------
O=[1;0]; l=[0;1];X=[0 1; 1 0]; I=eye(2); II=eye(2^2); III=eye(2^3);IIII=eye(2^4); IIIII=eye(2^5);
ll=kron(l,l)*kron(l',l'); OO = kron(O,O)*kron(O',O'); Ol=kron(O,l)*kron(O',l'); lO = kron(l,O)*kron(l',O'); 
ToffoliALL1 = kron(kron(kron(kron(kron(X,OO),I),I),I),I)+...
              kron(kron(kron(kron(kron(I,ll),I),I),X),I)+...
              kron(kron(kron(kron(kron(I,Ol),X),I),I),I)+...
              kron(kron(kron(kron(kron(I,lO),I),X),I),I);
ToffoliALL2 = kron(kron(kron(kron(kron(I,OO),I),I),X),I)+...
              kron(kron(kron(kron(kron(X,ll),I),I),I),I)+...
              kron(kron(kron(kron(kron(I,Ol),I),X),I),I)+...
              kron(kron(kron(kron(kron(I,lO),X),I),I),I);
% Coin = expm(-1i*pi/2*(kron(X,I)+kron(I,X))); Coin = kron(Coin,IIIII);

Coin = expm(-1i*pi/4*(kron(X,I)+kron(I,X))); Coin = kron(kron(I,Coin),IIII);

Par.Utarg=ToffoliALL2*Coin*ToffoliALL1*Coin;
Par.SaveFileName='Test';
Par.GuessFlag=1;
Par.GuessFileName='CR7_ToffoliALL_C2Step2d';
Par.nSec = 38;
Par.TargFid = 0.999;
% %---------------------------------------------------------------------
% O=[1;0]; l=[0;1];X=[0 1; 1 0]; I=eye(2); II=eye(2^2); III=eye(2^3);IIIII=eye(2^5);
% ll=kron(l,l)*kron(l',l'); OO = kron(O,O)*kron(O',O'); Ol=kron(O,l)*kron(O',l'); lO = kron(l,O)*kron(l',O'); 
% ToffoliALL1 = kron(kron(kron(kron(kron(X,OO),I),I),I),I)+...
              % kron(kron(kron(kron(kron(I,ll),I),I),X),I)+...
              % kron(kron(kron(kron(kron(I,Ol),X),I),I),I)+...
              % kron(kron(kron(kron(kron(I,lO),I),X),I),I);
% ToffoliALL2 = kron(kron(kron(kron(kron(I,OO),I),I),X),I)+...
              % kron(kron(kron(kron(kron(X,ll),I),I),I),I)+...
              % kron(kron(kron(kron(kron(I,Ol),I),X),I),I)+...
              % kron(kron(kron(kron(kron(I,lO),X),I),I),I);
% % Coin = expm(-1i*pi/10*(kron(X,I)+kron(I,X))); Coin = kron(Coin,IIIII);
% Par.Utarg=ToffoliALL2*ToffoliALL1;
% Par.SaveFileName='CR7_ToffoliALL_Step2a';
% Par.GuessFlag=1;
% Par.GuessFileName='CR7_ToffoliALL_Step1a';
% Par.nSec = 17;
% %---------------------------------------------------------------------
% O=[1;0]; l=[0;1];X=[0 1; 1 0]; I=eye(2); II=eye(2^2); III=eye(2^3);
% ll=kron(l,l)*kron(l',l'); OO = kron(O,O)*kron(O',O'); Ol=kron(O,l)*kron(O',l'); lO = kron(l,O)*kron(l',O'); 
% ToffoliALL1 = kron(kron(kron(kron(kron(X,OO),I),I),I),I)+...
%               kron(kron(kron(kron(kron(I,ll),I),I),X),I)+...
%               kron(kron(kron(kron(kron(I,Ol),X),I),I),I)+...
%               kron(kron(kron(kron(kron(I,lO),I),X),I),I);
% Par.Utarg=ToffoliALL1;
% Par.SaveFileName='CR7_ToffoliALLg';
% Par.GuessFlag=1;
% Par.GuessFileName='CR7_ToffoliALLf';
% Par.nSec = 29;
% %---------------------------------------------------------------------
% O=[1;0]; l=[0;1];I=eye(2); II=eye(2^2); III=eye(2^3); IIII=eye(2^4);
% A=O*O';  B=O*l'; C=l*O';  D=l*l';
% X=[0 1; 1 0];
% % 231 236 234 235
% % 123 236 234 235
% U1 = kron(kron(kron(X,A),A),IIII) + kron(kron(kron(kron(kron(I,A),D),II),X),I)...
%     +kron(kron(kron(kron(I,D),A),X),III)+...
%      kron(kron(kron(kron(kron(I,D),D),I),X),II);
% Par.Utarg=U1;
% Par.SaveFileName='CR7_RandomWalkStep1a';
% Par.GuessFlag=1;
% Par.GuessFileName='CR7_RandomWalkStep1';
% Par.nSec = 35;
%---------------------------------------------------------------------

%---------------------------------------------------------------------
% load('SevQubitsCrotonicC4','U2');
% Par.Utarg=U2;
% Par.SaveFileName='CR7_Steane_U2';
% Par.GuessFlag=0;
% Par.GuessFileName='CR7_Steane_U1';
% Par.nSec = 14;
%---------------------------------------------------------------------
%---------------------------------------------------------------------
% load('SevQubitsCrotonic','Err','E3','E2','E1');
% Par.Utarg=Err;
% Par.SaveFileName='UErr_Cr8';
%---------------------------------------------------------------------
%---------------------------------------------------------------------
% Par.Utarg=Cnot_IJaa(2,4);
% trace(Par.Utarg*Par.Utarg')
% Par.SaveFileName='SWAP24';
%---------------------------------------------------------------------
%---------------------------------------------------------------------
% [Ix,Iy,Iz,~,~,~,~] = prodopSparse([1/2 1/2],[4 3]);
% Par.Utarg = expm(-1i*pi/2*(Ix(:,:,2)))*(exp(-1i*pi*Iz(:,2).*Iz(:,1)).*(exp(-1i*pi*Iz(:,2).*Iz(:,6)).*(expm(-1i*pi/2*(Ix(:,:,2)+Iy(:,:,3)))*(exp(-1i*pi*Iz(:,3).*Iz(:,2)).*(exp(-1i*pi*Iz(:,3).*Iz(:,5)).*(expm(-1i*pi/2*(Ix(:,:,4)+Iy(:,:,3)))*(exp(-1i*pi*Iz(:,4).*Iz(:,3)).*(expm(-1i*pi/2*(Ix(:,:,7)+Iy(:,:,4)))*(exp(-1i*pi*Iz(:,4).*Iz(:,7)).*expm(-1i*pi/2*Iy(:,:,7)))))))))));
% A1=Par.Utarg*diag(Iz(:,7))*Par.Utarg';
% A2=diag(Iz(:,1))*diag(Iz(:,2))*diag(Iz(:,3).*Iz(:,4).*Iz(:,5).*Iz(:,6).*Iz(:,7));
% trace(A1*A2)/sqrt(trace(A1^2)*trace(A2^2))
% sdf
% Par.SaveFileName='CR7_PPS_Encoding';
%---------------------------------------------------------------------
% [Ix,Iy,Iz,~,~,~,~] = prodopSparse([1/2 1/2],[4 3]);
% Par.Utarg = expm(1i*pi/2*Iy(:,:,4))*(exp(-1i*pi*Iz(:,4).*Iz(:,7)).*(expm(1i*pi/2*(Iy(:,:,3)-Ix(:,:,4)))*(exp(-1i*pi*Iz(:,4).*Iz(:,3)).*(expm(1i*pi/2*(Iy(:,:,5)-Ix(:,:,3)))*(exp(-1i*pi*Iz(:,5).*Iz(:,3)).*(expm(1i*pi/2*(Iy(:,:,2)-Ix(:,:,5)))*(exp(-1i*pi*Iz(:,2).*Iz(:,3)).*(expm(1i*pi/2*(Iy(:,:,6)-Ix(:,:,2)))*(exp(-1i*pi*Iz(:,2).*Iz(:,6)).*(expm(1i*pi/2*(Iy(:,:,1)-Ix(:,:,6)))*(exp(-1i*pi*Iz(:,2).*Iz(:,1)).*expm(-1i*pi/2*(Ix(:,:,1))))))))))))));
% % Rin=zeros(2^7); Rin(128,1)=1; Rin(1,128)=1;
% % A1=Par.Utarg*Rin*Par.Utarg';
% % Rout=zeros(2^6);Rout(1,1)=1;
% % A2=kron(Rout,[0 1;1 0]);
% % trace(A1*A2)/sqrt(trace(A1^2)*trace(A2^2))
% Par.SaveFileName='CR7_PPS_Decoding';
% Par.GuessFlag=1;
% Par.GuessFileName='CR7_PPS_Decoding1';
% Par.nSec = 11;
%---------------------------------------------------------------------
% cnot = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0];
% Par.Utarg = kron(eye(2),kron(cnot,eye(2^6)));
% Par.SaveFileName='tt';
%---------------------------------------------------------------------

%---------------------------------------------------------------------
% [Ix,Iy,Iz,~,~,~,~] = prodopSparse([1/2 1/2],[4 5]);
% a=deg2rad(68.2);
% U1 = expm(1i*pi/2*(Ix(:,:,7)+Ix(:,:,8)+Ix(:,:,9)))*...
%     expm(-1i*pi/2*(cos(a)*(Ix(:,:,7)+Ix(:,:,8)+Ix(:,:,9))+sin(a)*(Iy(:,:,7)+Iy(:,:,8)+Iy(:,:,9))));
% 
% Uevo = diag(exp(-1i*pi*(Iz(:,4).*Iz(:,7)+Iz(:,4).*Iz(:,8)+Iz(:,4).*Iz(:,9))));
% U2=expm(-1i*pi/2*Iy(:,:,4));
% U3=expm(1i*pi/2*Iy(:,:,4));
% 
% a=deg2rad(116.6);
% U4 = expm(-1i*pi/2*(cos(a)*(Ix(:,:,7)+Ix(:,:,8)+Ix(:,:,9))+sin(a)*(Iy(:,:,7)+Iy(:,:,8)+Iy(:,:,9))))*...
%     expm(-1i*pi/2*(Ix(:,:,7)+Ix(:,:,8)+Ix(:,:,9)));
% 
% Par.Utarg=U4*Uevo*U3*Uevo*U2*Uevo*U1;
% Par.SaveFileName='USpinSelect2';
% clear Ix Iy Iz
%---------------------------------------------------------------------
% % swap = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
% % swap = kron(kron(eye(2),swap),eye(2));
% Par.SaveFileName='CnotCr';
% Par.Utarg = cnot;
% U=eye(2^3);
% U(7,7)=0; U(8,8)=0;U(7,8)=1; U(8,7)=1;
% Par.Utarg=kron(U,eye(2^4));
% Par.SaveFileName='test7Full';

% cnot = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0];
% cnot23 = kron(kron(eye(2),cnot),eye(2));
% Had = [1 1; 1 -1]/sqrt(2);              ZZ = [1 0; 0 -1];
% Had = kron(kron(eye(2),Had),eye(4));    ZZ = kron(kron(eye(2),ZZ),eye(4));
% Par.Utarg = ZZ*cnot23*Had*cnot23;
% Par.SaveFileName='UAnAs';

% ExtrasPath = ['/Users/hemantkatiyar/Google Drive/Documents/MATLAB/Compiler_strong' filesep 'Extras/UBAn_PURPLE' filesep];
% load([ExtrasPath 'UBAn10.mat']);
% Par.SaveFileName='Pu_UBAn10s';
% Par.Utarg = kron(Expression1,eye(4));
% Target fidelity of decomposition needed
Par.TargFid = 0.999;


% Par.FreeEvoMin = 20e-3;
% Par.Weight = 0.1;
