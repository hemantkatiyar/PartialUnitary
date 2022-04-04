% Author : Hemant Katiyar, 30-June-2020, (hkatiyar@uwaterloo.ca)

% Function to help decompose complex unitaries into single qubit gates
% and free evolutions
%
% Usage: Udecompose(ParFile)
% ParFile is the name (string) of input file stored in InputFiles folder

function Udecompose2(ParFile)

% Definining Some Paths
MolPath = [pwd filesep 'MoleculeFiles' filesep ];
InputPath = [pwd filesep 'InputFiles' filesep ];
SavePath = [pwd filesep 'SaveOutputs' filesep 'SaveOutputsDecom' filesep];
ExtrasPath = [pwd filesep 'Extras' filesep];
ToolsPath = [pwd filesep 'Tools'];
MainFilesPath = [pwd filesep 'MainPrograms' filesep 'MainProgramsDecom'];
addpath(ToolsPath);
addpath(MainFilesPath);

% Run the user define input file and Molecule file
run([InputPath ParFile]);
run([MolPath Par.MoleculeFile]);


% Number of spins
Mol.nSpin = sum(Mol.spinlist);

% Generate Paulis
% [Ix,Iy,Iz,~,~,~,D] = prodopSparse(1/2,Mol.nSpin);---
[Ix,Iy,Iz,~,~,~,D] = prodopSparse(Mol.spinNumbers,Mol.spinlist);

% Create Hint(free of chemical shift) and HintFull(includes chemical shift)
vtemp = zeros(size(Mol.v));
Mol.Hint = genHintWeak(Mol.spinlist,vtemp,Mol.J,D,Iz);
Mol.HintFull = genHint(Mol.spinlist,Mol.v,Mol.J,D,Ix,Iy,Iz);

% Hadamard gates acting on all spins
% Had = (1/sqrt(2))*[1 1; 1 -1];
% for j=1:Mol.nSpin-1
%     Had = kron(Had,(1/sqrt(2))*[1 1; 1 -1]);
% end
Iyy=zeros(D);
for j=1:Mol.nSpin
    Iyy = Iyy + Iy(:,:,j);
end
Had = expm(-1i*pi/2*Iyy);

% Initializing some variables
Par.VarPerSec = Mol.nSpin*2 + 1;
Par.TotVar = Par.VarPerSec*Par.nSec + Mol.nSpin*3;

% Either create initial guess or load a guess depending upon input from user
if Par.GuessFlag==0
    x = InitialGuess(Par.TotVar,Par.nSec,Par.VarPerSec,Par.DelayControl,Par.MaxIniDelay);
%     save guess.mat x
%     load guess1.mat
else
    load([SavePath Par.GuessFileName],'-mat','x');
%      x = x+(rand(size(x))-0.5);
     % x = [rand(1,3*15)*5 x];
%     x= x+(rand(size(x)))/2;
%       x= x+(rand(size(x))-0.5);
end

% Calculate total delay time for guess
DelayTime=0;
for j=1:Par.nSec
    DelayTime=DelayTime+ abs(x(j*Par.VarPerSec))*Par.DelayControl/pi;
end

% Save data in savepath
fileName=[SavePath Par.SaveFileName '.mat'];
save(fileName,'Par','Mol','x')

% % Function Definition
% CalcFid = @(x,y,z) abs(trace(x*y'))/2^z;


%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN PROGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
stepsize = Par.init_stepsize;

[Fid,Uf] = CalcFidDecom(x,Mol.nSpin,Par.nSec,Par.VarPerSec,Par.DelayControl,Iz,Had,Mol.Hint,Par.Utarg,D);

oldgrad = zeros(1,Par.TotVar);
olddirc = zeros(1,Par.TotVar);

counter=1;
fprintf('  Fidelity  |  StepSize  |  Max(ConjGrad)  |  Iteration  |  FreeEvo(ms)  |  Time(s) \n');
fprintf('  %6.6f    %6.5f        %6.5f        %5d            %5.4f        %5.2f \n',Fid,stepsize,max(olddirc),counter,DelayTime*1e+3,toc);
dspCheck = max((Par.TargFid-Fid)*0.2,1e-4)+Fid;
FixPhases(Par.SaveFileName)
save(fileName,'Par','Mol','x','Fid')
% Conjugate gradient reset flag, if its 1 we lose conjugacy
resetflag=0;
% guesscount=0;

while Fid<Par.TargFid
    
    Grad = CalcGradDecom(x,Mol.nSpin,Par.nSec,Par.VarPerSec,Par.DelayControl,Par.Utarg,Ix,Iy,Iz,Mol.Hint,Uf,D);
    [gooddirc,multiFac] = ConjugateGradDecom(resetflag,Mol.nSpin,counter,Grad,oldgrad,olddirc,Fid,x,stepsize,Iz,Had,Mol.Hint,Par.nSec,Par.VarPerSec,Par.DelayControl,Par.Utarg,D);
    
    resetflag=0;
    
    oldgrad = Grad;
    oldstepsize = stepsize;
    oldx = x;
    oldFid = Fid;
    olddirc = gooddirc;
    
    x = x+multiFac*stepsize*gooddirc;
    stepsize = sqrt(multiFac)*stepsize;
    [Fid,Uf] = CalcFidDecom(x,Mol.nSpin,Par.nSec,Par.VarPerSec,Par.DelayControl,Iz,Had,Mol.Hint,Par.Utarg,D);
    
    % We screwed up, either fit was not accurate or conjugate gradient are bad
    % Do a more fine linsesearch
    if Fid<oldFid
        resetflag=1;
        mulFineRange = linspace(0,4,201);
        FidTemp = zeros(1,length(mulFineRange));
        for j=1:length(mulFineRange)
            xtemp = oldx+mulFineRange(j)*oldstepsize*Grad;
            FidTemp(j)= CalcFidDecom(xtemp,Mol.nSpin,Par.nSec,Par.VarPerSec,Par.DelayControl,Iz,Had,Mol.Hint,Par.Utarg,D);
        end
        [~,indx] = max(FidTemp);
        multiFacTemp = mulFineRange(indx);
        multiFac = max(0.01,multiFacTemp);
        x = oldx+multiFac*stepsize*Grad;   
        stepsize=oldstepsize*sqrt(multiFac);
        [Fid,Uf] = CalcFidDecom(x,Mol.nSpin,Par.nSec,Par.VarPerSec,Par.DelayControl,Iz,Had,Mol.Hint,Par.Utarg,D);
    end
    
    counter=counter+1;
    
    % Displaying and saving current variables
    if (Fid>dspCheck || mod(counter,50)==0)
        save(fileName,'Par','Mol','x','Fid')
        FixPhases(Par.SaveFileName)
        DelayTime=0;
        for j=1:Par.nSec
            DelayTime=DelayTime+abs(x(j*Par.VarPerSec))*Par.DelayControl/pi;
        end
        fprintf('  %6.6f    %6.5f        %6.5f        %5d            %5.4f          %5.2f \n',Fid,stepsize,max(olddirc),counter,DelayTime*1e+3,toc);
        dspCheck = max((Par.TargFid-Fid)*0.2,1e-4)+Fid;
    end
    
    % Check if the gradients become too small or there is not much improvement
    % in the fidelity. Depending on KeepTrying, we either stop and print out
    % the reason for stopping or start with a new guess
    if max(gooddirc)<1e-5 || abs(Fid-oldFid)<1e-9
        if Par.KeepTrying==0 %If user asked not to try new guesses
            if max(gooddirc)<1e-5
                fprintf('Try new guess, small gradients \n');
%                 Uf{end}
            elseif abs(Fid-oldFid)<1e-9
                fprintf('Try new guess, small fidelity improvement \n');
%                 Uf{end}
            end
            return
        else %If user asked to keep on trying the guess
            fprintf('Trying a new guess \n');
            x = InitialGuess(Par.TotVar,Par.nSec,Par.VarPerSec,Par.DelayControl,Par.MaxIniDelay);
            DelayTime=0;
            for j=1:Par.nSec
                DelayTime=DelayTime+abs(x(j*Par.VarPerSec))*Par.DelayControl/pi;
            end
            stepsize=Par.init_stepsize;
            [Fid,Uf] = CalcFidDecom(x,Mol.nSpin,Par.nSec,Par.VarPerSec,Par.DelayControl,Iz,Had,Mol.Hint,Par.Utarg,D);

            oldgrad = zeros(1,Par.TotVar);
            olddirc = zeros(1,Par.TotVar);
            counter=1;
            fprintf('  Fidelity  |  StepSize  |  Max(ConjGrad)  |  Iteration  |  FreeEvo(ms)  |  Time(s) \n');
            fprintf('  %6.6f    %6.5f        %6.5f        %5d            %5.4f        %5.2f \n',Fid,stepsize,max(olddirc),counter,DelayTime*1000,toc);
            dspCheck = max((Par.TargFid-Fid)*0.2,1e-4)+Fid;
        end
    end
end %Big while loop

save(fileName,'Par','Mol','x','Fid')
% Uf{end}
FixPhases(Par.SaveFileName)

%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN PROGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add penalty on time
% write correcting phases program auto ---- Done
% test with other stepsize function ---- Better results
% test new conjugate method ----- Bad results
% Add overwrite Prevention
