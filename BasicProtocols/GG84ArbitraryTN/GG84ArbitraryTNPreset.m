function qkdInput = GG84ArbitraryTNPreset()
% BasicBB84PersonalPreset
% EB-BB84 con canale ignoto (vincoli su rho_AB).
% - Description: G=Id, pinching Z su A direttamente su AB.
% - Channel: vincoli per rho_A = I/2 e QBER^X = QBER^Z = p.
% - Funzionale: D(rho || Z(rho)) (BasicKeyRateFunc).
% Salva anche parametri fisici e lista p per uso nel main.

qkdInput = QKDSolverInput();

%% Parameters
qkdInput.addFixedParameter("renyi", false);
qkdInput.addFixedParameter("alpha", 1.00000000001);
qkdInput.addFixedParameter("distance", 20);
qkdInput.addFixedParameter("eps", 1e-8);
qkdInput.addScanParameter("EveDisturbance", num2cell(linspace(0.01,0.15,12)));
qkdInput.addScanParameter("flipProb", num2cell([0, 0.2, 0.5]))
%qkdInput.addScanParameter("flipProb", num2cell(linspace(0.2,0.3,15)));
%qkdInput.addFixedParameter("flipProb", 0.3);

%% Modules
descriptionModule = QKDDescriptionModule(@GG84ArbitraryTNDescriptionFunc);
qkdInput.setDescriptionModule(descriptionModule);

channelModule = QKDChannelModule(@GG84ArbitraryTNChannelFunc);
qkdInput.setChannelModule(channelModule);

keyRateModule = QKDKeyRateModule(@GG84ArbitraryTNKeyRateFunc);
qkdInput.setKeyRateModule(keyRateModule);

%% Math solver
mathSolverOptions = struct();
mathSolverOptions.initMethod = 1;
mathSolverOptions.frankWolfeMethod = @FrankWolfe.vanilla;
mathSolverOptions.frankWolfeOptions = struct("maxIter",20,"maxGap",1e-5);
mathSolverOptions.linearConstraintTolerance = 1e-10;
mathSolverMod = QKDMathSolverModule(@GG84ArbitraryTNFW2StepSolver, mathSolverOptions);
qkdInput.setMathSolverModule(mathSolverMod);

%% Opzioni globali
qkdInput.setGlobalOptions(struct("errorHandling", ErrorHandling.CatchWarn, ...
                                 "verboseLevel", 1, ...
                                 "cvxSolver", "mosek"));
end