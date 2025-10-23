function qkdInput = BasicBB84Alice2DFinitePresetOpt(alpha)
% BasicBB84Alice2DFinitePreset A preset for BB84 with qubits, for
% finite-size keyrates based upon Phys. Rev. Research 3, 013274

qkdInput = QKDSolverInput();

%% Parameters

%Add loss
lossdB = linspace(0,0,1);
transmittance = 10.^(-lossdB/10);
qkdInput.addScanParameter("transmittance", num2cell(transmittance));

% qkdInput.addScanParameter("numSignals",num2cell(10.^linspace(5,8,4)));
%qkdInput.addScanParameter("misalignmentAngle", num2cell(pi/180*linspace(0,5,6)) );
qkdInput.addFixedParameter("pz",1/2);
qkdInput.addFixedParameter("fEC",1.2); %efficiency of error-correction.

qkdInput.addFixedParameter("depolarization",0.01);
% qkdInput.addFixedParameter("transmittance",0);

qkdInput.addFixedParameter("alpha", alpha);


%finite size parameters.

qkdInput.addFixedParameter("alphabetSize", 2); % encoding alphabet size; for qubits, this is 2

% epsilon

qkdInput.addFixedParameter("epsilonPE",(1/4)*1e-8); % acceptance test (parameter estimation) epsilon
qkdInput.addFixedParameter("epsilonPA",(1/4)*1e-8); % PA epsilon
qkdInput.addFixedParameter("epsilonEC",(1/4)*1e-8); % error-verification epsilon
qkdInput.addFixedParameter("epsilonBar",(1/4)*1e-8);% smoothing epsilon



qkdInput.addFixedParameter("numSignals",3e5); %total number of signals
qkdInput.addFixedParameter("misalignmentAngle", 0); %pi*8.1301/180
qkdInput.addFixedParameter("pTest", 0.2); %m=pTest*numSignals used for testing






qkdInput.addFixedParameter("tExp", -7);




% description is the same as the asymptotic qubit BB84
% descriptionModule = QKDDescriptionModule(@BasicBB84Alice2DDescriptionFunc);
descriptionModule = QKDDescriptionModule(@BasicBB84LossyDescriptionFunc);
qkdInput.setDescriptionModule(descriptionModule);

% channel model 

% channelModule = QKDChannelModule(@BasicBB84Alice2DChannelFunc);
channelModule = QKDChannelModule(@BasicBB84LossyChannelFunc);
qkdInput.setChannelModule(channelModule);
    
% key rate module
keyRateOptions = struct();
keyMod = QKDKeyRateModule(@BasicBB84Alice2DFiniteKeyRateFunc, keyRateOptions);
qkdInput.setKeyRateModule(keyMod);

optimizerMod = QKDOptimizerModule(@coordinateDescentFunc,struct("verboseLevel",0),struct("verboseLevel",0));
qkdInput.setOptimizerModule(optimizerMod);

% math solver options
mathSolverOptions = struct();
mathSolverOptions.initMethod = 1;
mathSolverOptions.maxIter = 10;
mathSolverOptions.maxGap = 1e-6;
mathSolverOptions.blockDiagonal = false;
mathSolverMod = QKDMathSolverModule(@FW2StepSolver,mathSolverOptions,mathSolverOptions);
qkdInput.setMathSolverModule(mathSolverMod);

% global options
qkdInput.setGlobalOptions(struct("errorHandling",3,"verboseLevel",1,"cvxSolver","mosek","cvxPrecision","high"));

