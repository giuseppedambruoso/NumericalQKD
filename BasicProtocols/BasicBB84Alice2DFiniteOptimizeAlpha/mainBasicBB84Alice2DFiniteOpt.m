function [kr] = mainBasicBB84Alice2DFiniteOpt(alpha)

%This function is to obtain a more precise optimal alpha. To obtain a
%rough value, one could also use the "scanparameter" in preset
%First change the block size in file BasicBB84Alice2DFinitePresetOpt
%Then run: 
% >> ops = optimset('MaxIter',15, 'MaxFunEvals',15)
% >> [a_opt,kr]=fmincon(@(a) -mainBasicBB84Alice2DFiniteOpt(a),1.01,[],[],[],[],1.0001,a_opt,[],ops)


%% Begin here:
    %pick the preset file 
    disp(alpha)
    qkdInput = BasicBB84Alice2DFinitePresetOpt(alpha);

    %run the QKDSolver with this input
    results = MainIteration(qkdInput);

    %results is a cell containing different N values
    kr = results(1).keyRate;
    
    
    %save the results and preset to a file.
    
    save("BasicBB84Alice2DFiniteresultsOpt.mat","results","qkdInput");
    
    
    %% plot the result
    %QKDPlot.simple1DPlot(qkdInput,results)
end 