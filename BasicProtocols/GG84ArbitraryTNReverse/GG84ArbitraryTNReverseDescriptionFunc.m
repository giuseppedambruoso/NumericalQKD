function [newParams,modParser] = GG84ArbitraryTNReverseDescriptionFunc(params, options,debugInfo)
% BasicBB84Alice2DDescriptionFunc A simple description function for a qubit
% BB84 protocol with no loss that uses the Schmidt decomposition to turn
% Here, Schmidt decomposition was used to shrink Alice from a 4d space to a
% 2d space. This is used with BasicKeyRateFunc.
%
% Input parameters:
% * pz: The probability that Alice measures in the Z-basis (for this
%   protocol, it's also the probability that Bob measures in the Z-basis as
%   well). It must be between 0 and 1.
% Output parameters:
% * observablesJoint: The joint observables for Alice and Bob's measurement
%   of the signals.
% * dimA: dimension of Alice's system.
% * dimB: dimension of Bob's system.
% * rhoA: Alice's reduced density matrix for prepare and measure based
%   protocols.
% * POVMA: Alice's set of POVM operators which she measures her state with
%   in the source replacement scheme.
% * POVMB: Bob's set of POVM operators which he measures his state with.
% * announcementsA: Alice's announcements for each of her POVM operators.
%   Can be integers or strings.
% * announcementsB: Bob's announcements for each of his POVM operators.
%   Can be integers or strings.
% * keyMap: An array of KeyMap objects that contain pairs of accepted
%   announcements and an array dictating the mapping of Alice's measurement
%   outcome to key bits (May be written with Strings).
% * krausOps: A cell array of matrices. The Kraus operators that form the G
%   map on Alice and Bob's joint system. These should form a completely
%   positive trace non-increasing linear map. Each Kraus operator must be
%   the same size.
% * keyProj:  A cell array of projection operators that extract the key
%   from G(\rho). These projection operators should sum to identity. This
%   map is often called Z.
% Options:
% * none
% DebugInfo:
% * krausSum: sum_i K^\dagger_i*K_i which should be <= identity for
%   a CPTNI map.
%
% See also QKDDescriptionModule, BasicKeyRateFunc
arguments
    params (1,1) struct
    options (1,1) struct
    debugInfo (1,1) DebugInfo
end

%% options parser
optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.parse(options);
options = optionsParser.Results;
%% module parser
modParser = moduleParser(mfilename);
modParser.addRequiredParam("flipProb",...
    @isscalar,...
    @(x) mustBeInRange(x,0,1))
modParser.parse(params)
params = modParser.Results;

%% simple setup
newParams = struct();
q = params.flipProb;

%% dimension sizes of Alice and Bob
dimA = 2;
dimB = 2;

%% rhoA for source replacement scheme
newParams.rhoA = eye(dimA)/dimA;

%% Kraus Ops (for G map)
% A: Alice's system, B: Bob's System, C: announcement register, R:
% key register.
% The Kraus operators are matrices from ABC \rightarrow RBC. Here we used
% an isometry to shrink the Kraus operators from outputting on RABC to just
% RBC. This lets us save time on computing eigen values later. The factor
% of pz comes from a \sqrt(pz) from Alice's measurements(from Schmidt

% Reverse reconciliation
ketZ = [1;0]; ketO = [0;1];
krausOps = {kron(eye(dimA),kron(sqrt(1-q)*eye(dimB)+sqrt(q)*(ketZ*ketO'+ketO*ketZ'),zket(2,1)))};
%Here we compute sum_i K^\dagger_i*K_i. Which should satisfy sum_i
%K^\dagger_i*K_i <= I. A.K.A. the Kraus operators represent a
%completely positive, trace non-increasing linear map.
krausSum = 0;
for index = 1:numel(krausOps)
    krausSum = krausSum+krausOps{index}'*krausOps{index};
end
debugInfo.storeInfo("krausSum",krausSum);

%% key projection
% the key projection is formed from projection operators on the space RBC. The
% projectors form a pinching map that measures the key register.
% (in a sense, we want to break any quantum correlations between this
% register and any of the remaining registers.)
proj0 = kron(eye(dimB),kron(diag([1,0]),eye(dimB)));
proj1 = kron(eye(dimB),kron(diag([0,1]),eye(dimB)));
keyProj = {proj0,proj1};

%% set Kraus ops and key projection in new parameters
newParams.krausOps = krausOps;
newParams.keyProj = keyProj;

end