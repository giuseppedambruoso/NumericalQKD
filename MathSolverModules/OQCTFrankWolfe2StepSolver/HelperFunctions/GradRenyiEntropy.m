function gradobj = GradRenyiEntropy(alpha,gRho,zRho, krausOperators,keyProj)
%First, get eigenvectors and eigenvalues of gRho
% For consistency get the same perturbation value for both.
[PG,DG] = eig(gRho);

%Get eigenvectors and eigenvalues of zRho
[PZ,DZ] = eig(zRho);

Sz = size(gRho);
grad1 = zeros(Sz);
grad2 = zeros(Sz);
grad3 = zeros(Sz);
beta = 1/alpha;
mu = (1-beta)/(2*beta);
Xi = (zRho^mu) * gRho * zRho^mu;

for k = 1:1:Sz(1)
    for l = 1:1:Sz(1)
         grad1 = grad1+gradfrac(DZ(k,k),DZ(l,l), mu)* PZ(:,k)*PZ(:,k)'...
             *gRho*(zRho^mu)*(beta*Xi^(beta-1))*PZ(:,l)*PZ(:,l)';
         grad3 = grad3+gradfrac(DZ(k,k),DZ(l,l),mu)*PZ(:,k)*PZ(:,k)'...
             *(beta*Xi^(beta-1))*(zRho^mu)*gRho*PZ(:,l)*PZ(:,l)';
    end
end
grad1 = ApplyMap(grad1,keyProj);
grad2 = (zRho^mu)*(beta*Xi^(beta-1))*zRho^mu;
grad3 = ApplyMap(grad3,keyProj);


gradQ = grad1+grad2+grad3;

%% Sandwiched divergence
Q = trace(Xi^beta);
gradD = (1/(beta-1))*ApplyMap((gradQ)/Q-eye(Sz)/trace(gRho),DualMap(krausOperators));
gradobj = ApplyMap(eye(Sz)*RenyiEntropy(alpha,gRho,zRho)/trace(gRho)...
    , DualMap(krausOperators))+trace(gRho)*gradD;




end

function gradient = gradfrac(ek,el,power)
    if (ek == el || abs(ek-el)<= 1e-10)
    % if abs(ek-el)<1e-10
        gradient=power*ek^(power-1);
    else
        gradient = (ek^power - el^power)/(ek-el);
    end
end
% =======
