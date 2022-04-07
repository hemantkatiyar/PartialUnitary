;Par=[];
% Name of molecule file located in MoleculeFiles folder
Par.MoleculeFile='MolCrotonic';

% Number of sections in decomposition
Par.nSec = 0;

% Max value of sum of all delay times for the guess (in ms)
Par.MaxIniDelay = 100;

% Parmeter to control t
Par.DelayControl = 0.001;
% Initial step size(usually left alone)
Par.init_stepsize = 1;

% If GuessFlag=0, then nohe delay times(usually left alone)
% Par. guess is used, if its 1 load the GuessFileName
Par.GuessFlag=0;
Par.GuessFileName='ss';
% If algorithm stops prematurely, should we keep trying? if yes, then 1
Par.KeepTrying=0;

% Save the results wiht the following name in SaveOutputs\SaveOutputsDecom
Par.SaveFileName='Test';

% Define Unitary that needs decomposing


[Ix,Iy,Iz,~,~,~,~] = prodopSparse([1/2],[4]);
A=expm(-1i*pi/2*Ix(:,:,1));
Par.Utarg=A(:,1:2);
Par.TargFid = 0.999;
