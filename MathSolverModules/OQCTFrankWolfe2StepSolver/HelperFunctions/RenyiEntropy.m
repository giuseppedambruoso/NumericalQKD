function renyi_ent = RenyiEntropy(alpha,gRho,zRho)

beta = 1/alpha;
mu = (1-beta)/(2*beta);
Trace = trace(((zRho^mu)*gRho*zRho^mu)^beta);

renyi_ent = trace(gRho)*(1/(beta-1))*(log(Trace/trace(gRho)));
end

