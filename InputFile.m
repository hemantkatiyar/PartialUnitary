;Par=[];
% Name of molecule file located in MoleculeFiles folder
Par.MoleculeFile='MolCrotonic7';

% Number of sections in decomposition
Par.nSec = 10;

% Max value of sum of all delay times for the guess (in ms)
Par.MaxIniDelay = 100;

% Parmeter to control t
Par.DelayControl = 0.001;
% Initial step size(usually left alone)
Par.init_stepsize = 1;

% If GuessFlag=0, then nohe delay times(usually left alone)
% Par. guess is used, if its 1 load the GuessFileName
Par.GuessFlag=0;
Par.GuessFileName='CR7_PPS_Encoding10';
% If algorithm stops prematurely, should we keep trying? if yes, then 1
Par.KeepTrying=0;

% Save the results wiht the following name in SaveOutputs\SaveOutputsDecom
Par.SaveFileName='Test';

% Define Unitary that needs decomposing


% [Ix,Iy,Iz,~,~,~,~] = prodopSparse([1/2],[4]);
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

Coin = expm(-1i*pi/4*(kron(X,I)+kron(I,X))); Coin = kron(kron(I,Coin),IIII);

Par.Utarg=ToffoliALL2*Coin*ToffoliALL1*Coin;
Par.TargFid = 0.999;
