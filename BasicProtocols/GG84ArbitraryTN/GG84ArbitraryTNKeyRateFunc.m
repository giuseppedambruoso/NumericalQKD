function [keyRate, modParser] = GG84ArbitraryTNKeyRateFunc(params,options,mathSolverFunc,debugInfo)
% BasicKeyRateFunc A simple key rate function for a asymptotic key rate
% calculations with equality constraints and deterministic key map.
%
% Input parameters:
% * dimA: dimension of Alice's system.
% * dimB: dimension of Bob's system.
% * announcementsA: Array of announcements made for each measurement Alice
%   made ordered in the same way as the columns of expectationsJoint.
% * announcementsB: Array of announcements made for each measurement Bob
%   made ordered in the same way as the rows of expectationsJoint.
% * keyMap: An array of KeyMapElement objects that contain pairs of
%   accepted announcements and an array dictating the mapping of Alice's
%   measurement outcome to key bits (May be written with Strings).
% * krausOps: A cell array of matrices. The Kraus operators that form the G
%   map on Alice and Bob's joint system. These should form a completely
%   positive trace non-increasing linear map. Each Kraus operator must be
%   the same size.
% * keyProj:  A cell array of projection operators that perform the
%   pinching map key on  G(\rho). These projection operators should sum to
%   identity.
% * fEC: error correction efficiency. Set to 1 means for Shannon limit.
% * observablesJoint: The joint observables of Alice and Bob's
%   measurements. The observables must be Hermitian and each must be the
%   size dimA*dimB by dimA*dimB. The observables assume the spaces are
%   ordered A \otimes B. They also should be positive semi-definite and
%   should sum to identity.
% * expectationsJoint: The joint expectations (as an array) from Alice and
%   Bob's measurements that line up with it's corresponding observable in
%   observablesJoint. These values should be between 0 and 1.
% * rhoA (nan): The fixed known density matrix on Alice's side for
%   prepare-and-measure protocols. Setting to nan means this is ignored and
%   no rhoA constraint is passed to the convex solver. Furthermore, an
%   explicit Tr[rhoAB] == 1 constraint is added in its place.
% Outputs:
% * keyrate: Key rate of the QKD protocol measured in bits per block
%   processed.
% Options:
% * verboseLevel: (global option) See makeGlobalOptionsParser for details.
% DebugInfo:
% * deltaLeak: Error correction cost calculated for the protocol.
% * keyRateRelEntStep2Linearization: Estimation of the key rate by using
%   the relative entropy at the point where the Frank-Wolfe solver starts
%   its step 2 linearization. THIS IS NOT A SAFE LOWER BOUND.
%
% Reviewed by Devashish Tupkary 2023/09/18
% See also QKDKeyRateModule, makeGlobalOptionsParser
arguments
    params (1,1) struct
    options (1,1) struct
    mathSolverFunc (1,1) function_handle
    debugInfo (1,1) DebugInfo
end

%% options parser
optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.parse(options);
options = optionsParser.Results;
%% modParser
modParser = moduleParser(mfilename);

modParser.addRequiredParam("renyi", @mustBeBoolean);
modParser.addRequiredParam("alpha", @(x) mustBeGreaterThan(x,1));
modParser.addRequiredParam("krausOps"); %, @isCPTNIKrausOps);
modParser.addRequiredParam("keyProj", @(x) mustBeAKeyProj(x));
modParser.addRequiredParam("EveDisturbance", @isscalar, @(x) mustBeInRange(x,0,1));
modParser.addOptionalParam("rhoA", nan, @(x) isequaln(x,nan) || isDensityOperator(x));
modParser.addRequiredParam("observablesJoint",@(x) allCells(x,@ishermitian));
modParser.addRequiredParam("expectationsJoint",@(x) mustBeInRange(x,-1,1));
modParser.addRequiredParam("flipProb", @(x) mustBeInRange(x,0,1));
modParser.addAdditionalConstraint(@isEqualSize,["observablesJoint","expectationsJoint"]);

modParser.parse(params);


params = modParser.Results;

%% simple setup
debugMathSolver = debugInfo.addLeaves("mathSolver");
mathSolverInput = struct();

%% Error correction
qber = params.EveDisturbance+params.flipProb-2*params.EveDisturbance*params.flipProb;
deltaLeak = -qber*log2(qber)-(1-qber)*log2(1-qber);
debugInfo.storeInfo("deltaLeak",deltaLeak);


%% translate for the math solver
%now we have all the parts we need to get a key rate from the a math
%solver, but we need to put it into a form it can understand.
%first we give it the Kraus operators for the G map and the projection
%operators for the key map (Z).
mathSolverInput.krausOps = params.krausOps;
mathSolverInput.keyProj = params.keyProj;
mathSolverInput.alpha = params.alpha;
mathSolverInput.renyi = params.renyi;

numObs = numel(params.observablesJoint);

mathSolverInput.equalityConstraints = arrayfun(@(x)...
    EqualityConstraint(params.observablesJoint{x},params.expectationsJoint(x)),1:numObs);

% Include rhoA (if given) or add an explicit Tr[rhoAB] == 1 constraint.
if ~isequaln(params.rhoA,nan)
    mathSolverInput.rhoA = params.rhoA;
else
    mathSolverInput.equalityConstraints(end+1) = EqualityConstraint(eye(params.dimA*params.dimB),1);
end

% now we call the math solver function on the formulated inputs, with
% it's options.
[relEnt,~] = mathSolverFunc(mathSolverInput,debugMathSolver);

%store the key rate
keyRate = max(0,relEnt-deltaLeak);

if options.verboseLevel>=1
    %ensure that we cut off at 0 when we display this for the user.
    fprintf("Key rate: %e\n",max(keyRate,0));
end

%set the linearization estimate key rate as well for debugging
if isfield(debugMathSolver.info,"relEntStep2Linearization")
    keyRateStep2Linearization = debugMathSolver.info.relEntStep2Linearization - deltaLeak;
    debugInfo.storeInfo("keyRateRelEntStep2Linearization",keyRateStep2Linearization)

    if options.verboseLevel>=2
        fprintf("Key rate using step 2 linearization intial value: %e\n",max(keyRateStep2Linearization,0))
    end
end
end