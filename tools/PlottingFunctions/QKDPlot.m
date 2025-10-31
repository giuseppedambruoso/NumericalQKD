classdef QKDPlot
    % QKDPlot Helper class for plotting QKD results.
    %
    % Supports 1D plots (single scan parameter vs key rate) and
    % 2D plots (key rate vs EveDisturbance for different flipProb values).
    properties (Constant)
        keyRateTag = "_keyRate_"; % tag for key rate plotting
    end

    methods(Static)
%% --- 1D Plot ---
function simple1DPlot(qkdInput,results,options)
    arguments
        qkdInput (1,1) QKDSolverInput{QKDSolverInputMustHaveExactlyOneScanParam}
        results (:,1) struct
        options.xScaleStyle (1,1) string {mustBeMember(options.xScaleStyle,["linear","log","dB"])} = "linear";
        options.yScaleStyle (1,1) string {mustBeMember(options.yScaleStyle,["linear","log","dB"])} = "linear";
        options.markerLineStyles (1,1) string = "x-";
        options.figAxis (1,1) matlab.graphics.axis.Axes = axes(figure());
    end

    % --- PLOT dei risultati numerici ---
    scanParamName = string(fieldnames(qkdInput.scanParameters));
    QKDPlot.plotParameters({results}, scanParamName, ...
        QKDPlot.keyRateTag, ...
        "xScaleStyle", options.xScaleStyle, ...
        "yScaleStyle", options.yScaleStyle, ...
        "markerLineStyles", options.markerLineStyles, ...
        "figAxis", options.figAxis, ...
        "addLegend", false);

    hold(options.figAxis, "on");

    % --- Aggiunge la curva teorica 1 - h2(p) ---
    h2 = @(p) (p==0 | p==1).*0 + ...
              (p>0 & p<1).*(-p.*log2(p) - (1-p).*log2(1-p));

    % Range x teorico: adatta in base ai dati reali
    currentParams = [results.currentParams];
    paramValues = [currentParams.(scanParamName)];
    p_th = linspace(min(paramValues), max(paramValues), 500);

    % Calcola curva teorica
    r_th = max(0,1 - 2*h2(p_th));

    % Disegna curva teorica
    plot(options.figAxis, p_th, r_th, 'k--', 'LineWidth', 1.8, ...
         'DisplayName', '1 - h_2(p) (theoretical)');

    hold(options.figAxis, "off");

    % --- Etichette e legenda ---
    xlabel(options.figAxis, scanParamName);
    ylabel(options.figAxis, 'Key Rate');
    legend(options.figAxis, 'show', 'Location', 'best');
    grid(options.figAxis, 'on');
    box(options.figAxis, 'on');
end
%% --- 2D Plot (key rate vs EveDisturbance for flipProb) ---
function simple2DPlot(qkdInput, results, options)
    arguments
        qkdInput (1,1)
        results (:,1) struct
        options.xScaleStyle (1,1) string {mustBeMember(options.xScaleStyle, ["linear","log","dB"])} = "linear"
        options.yScaleStyle (1,1) string {mustBeMember(options.yScaleStyle, ["linear","log","dB"])} = "log"
        options.markerLineStyles (:,1) string = ["x-", "o-", "*-", "s-", "d-"]
        options.figAxis (1,1) matlab.graphics.axis.Axes = axes(figure())
    end

    hold(options.figAxis, 'on');

    % Extract keyRate, EveDisturbance, flipProb
    keyRates = [results.keyRate];
    eveDists = [results.currentParams];
    eveDists = [eveDists.EveDisturbance];
    flipProbs = [results.currentParams];
    flipProbs = [flipProbs.flipProb];

    % Unique flipProb values
    uniqueFlipProbs = unique(flipProbs);

    % Entropia binaria (robusta ai bordi p=0 o p=1)
    h2 = @(p) ...
        (p==0 | p==1).*0 + ...
        (p>0 & p<1).*(-p.*log2(p) - (1-p).*log2(1-p));

    % Parametri Q e q
    Q = linspace(0, 0.11, 1000);
    qs = [0.00000001, 0.2, 0.5];

    % Funzione s(q,Q) con protezione contro piccoli negativi numerici
    s_fun = @(q, Qv) sqrt( max(0, 1 - 16*q*(1-q).*Qv.*(1-Qv)) );

    % r_infty(Q,q)
    rInf = @(Qv, q) ( 1 - h2(Qv) + h2( (1 + s_fun(q, Qv))/2 ) ) ...
                    - h2( Qv + q - 2*q.*Qv );

    % --- Plotta tutte le curve teoriche (stessa figura)
    for q = qs
        plot(options.figAxis, Q, rInf(Q, q), 'LineWidth', 1.8, ...
             'DisplayName', sprintf('r_\\infty(q=%.2g)', q));
    end

    % --- Plotta i risultati numerici per ogni flipProb
    for i = 1:length(uniqueFlipProbs)
        q = uniqueFlipProbs(i);

        % Select indices for this flipProb
        idx = flipProbs == q;

        % Select data for this flipProb
        xData = eveDists(idx);
        yData = keyRates(idx);

        % Plot curve
        plot(options.figAxis, xData, yData, ...
            options.markerLineStyles(mod(i-1,length(options.markerLineStyles))+1), ...
            'LineWidth', 1.2, ...
            'DisplayName', ['flipProb = ' num2str(q)]);
    end

    hold(options.figAxis, 'off');
    xlabel(options.figAxis, 'EveDisturbance');
    ylabel(options.figAxis, 'Key Rate');
    legend(options.figAxis, 'show', 'Location', 'best');
    grid(options.figAxis, 'on');
    box(options.figAxis, 'on');
end
        %% --- Plot parameters from files ---
        function plotParametersFromFiles(fileNames,xParamName,yParamName,options)
            arguments
                fileNames (:,1) string{mustBeFile}
                xParamName (1,1) string
                yParamName (1,1) string
                options.addLegend (1,1) logical =true;
                options.legendNames (:,1) string {mustBeEqualSize(options.legendNames,fileNames)} = fileNames;
                options.xScaleStyle (1,1) string {mustBeMember(options.xScaleStyle,["linear","log","dB"])} = "linear";
                options.yScaleStyle (1,1) string {mustBeMember(options.yScaleStyle,["linear","log","dB"])} = "linear";
                options.markerLineStyles (:,1) string = ["x-", "o-", ".-", "^-", "+-", ">-", "s-", "p-", "*-"];
                options.figAxis (1,1) matlab.graphics.axis.Axes = axes(figure());
            end

            resultsSets = cell(size(fileNames));
            for index = 1:numel(fileNames)
                tempFile = load(fileNames(index));
                resultsSets{index} = tempFile.results;
            end

            QKDPlot.plotParameters(resultsSets,xParamName,yParamName,...
                "legendNames",options.legendNames,"xScaleStyle",options.xScaleStyle,...
                "yScaleStyle",options.yScaleStyle,"markerLineStyles",options.markerLineStyles,...
                "figAxis",options.figAxis, "addLegend",options.addLegend);
        end

        %% --- Generic plotParameters ---
        function plotParameters(resultsSets,xParamName,yParamName,options)
            arguments
                resultsSets (:,1) cell {mustBeCellOf(resultsSets,"struct")}
                xParamName (1,1) string {mustBeANamedParameterOrKeyFlag(xParamName,resultsSets)}
                yParamName (1,1) string {mustBeANamedParameterOrKeyFlag(yParamName,resultsSets)}
                options.addLegend (1,1) logical =true;
                options.legendNames (:,1) string {mustBeEqualSize(options.legendNames,resultsSets)} = compose("data set %d",1:numel(resultsSets));
                options.xScaleStyle (1,1) string {mustBeMember(options.xScaleStyle,["linear","log","dB"])} = "linear";
                options.yScaleStyle (1,1) string {mustBeMember(options.yScaleStyle,["linear","log","dB"])} = "linear";
                options.markerLineStyles (:,1) string = ["x-", "o-", ".-", "^-", "+-", ">-", "s-", "p-", "*-"];
                options.figAxis (1,1) matlab.graphics.axis.Axes = axes(figure());
            end

            xUseDB = options.xScaleStyle == "dB";
            yUseDB = options.yScaleStyle == "dB";

            [xDataSets,simpleIndexingXFlag] = extractAndFormatDataSets(xParamName,resultsSets);
            [yDataSets,simpleIndexingYFlag] = extractAndFormatDataSets(yParamName,resultsSets);

            if simpleIndexingXFlag && xUseDB || simpleIndexingYFlag && yUseDB
                throw(MException("plotParameters:CantSimpleIndexAndDBScale",...
                    "The data requires simple linear indexing which is not compatible with using dB scaling style."))
            end
            if xUseDB
                xDataSets = cellfun(@(x)-10*log10(x),xDataSets,"UniformOutput",false);
            end
            if yUseDB
                yDataSets = cellfun(@(x)-10*log10(x),yDataSets,"UniformOutput",false);
            end

            xLabelName = formatLabelName(xParamName,simpleIndexingXFlag,xUseDB);
            yLabelName = formatLabelName(yParamName,simpleIndexingYFlag,yUseDB);

            plotDataSets(options.figAxis,xDataSets,yDataSets,xLabelName,yLabelName,...
                options.xScaleStyle~="log",options.yScaleStyle~="log",...
                options.markerLineStyles,options.legendNames,options.addLegend);
        end
    end
end

%% --- Helper functions (outside class) ---
function [dataSets, simpleIndexingFlag] = extractAndFormatDataSets(paramName,resultsSets)
dataSets = cell(size(resultsSets));
simpleIndexingFlag = false;

for index =1:numel(resultsSets)
    dataSets{index} = extractParameter(paramName,resultsSets{index});
    simpleIndexingFlag = ~all(cellfun(@(x)isscalar(x) && isreal(x), dataSets{index}));
    if simpleIndexingFlag, break; else, dataSets{index} = cell2mat(dataSets{index}); end
end

if simpleIndexingFlag
    numelEach = cellfun(@numel,resultsSets);
    if all(numelEach == numelEach(1))
        dataSets{:} = 1:numelEach(1);
    else
        throw(MException("plotParameters:SimleAssignmentOnlyForSameSize",...
            "If the parameters can't be converted to the single real values, then each result must be the same size."))
    end
end
end

function labelName = formatLabelName(paramName,simpleIndexingFlag,useDB)
if paramName == QKDPlot.keyRateTag, paramName = "key rate"; end
paramName = replace(paramName,"_","\_");
if simpleIndexingFlag
    labelName = paramName+" index";
elseif useDB
    labelName = paramName+" parameterized in dB";
else
    labelName = paramName;
end
end

function paramValues = extractParameter(paramName,results)
if paramName == QKDPlot.keyRateTag
    paramValues = {results(:).keyRate};
else
    currentParams = [results(:).currentParams];
    paramValues = {currentParams.(paramName)};
end
end

function plotDataSets(figAxis,xDataSets,yDataSets,xLabelName,yLabelName,linX,linY,markerLineStyles,legendNames, addLegend)
arguments
    figAxis (1,1) matlab.graphics.axis.Axes
    xDataSets (:,1) cell {mustBeCellOf(xDataSets,"double")}
    yDataSets (:,1) cell {mustBeCellOf(yDataSets,"double"),mustBeEqualSize(yDataSets,xDataSets)}
    xLabelName (1,1) string
    yLabelName (1,1) string
    linX (1,1) logical
    linY (1,1) logical
    markerLineStyles (:,1) string
    legendNames (:,1) string
    addLegend (1,1) logical
end

if linX, figAxis.XScale = 'linear'; else, figAxis.XScale = 'log'; end
if linY, figAxis.YScale = 'linear'; else, figAxis.YScale = 'log'; end

figAxis.XLabel.String = xLabelName;
figAxis.YLabel.String = yLabelName;

for index = 1:numel(xDataSets)
    if index == 1, hold(figAxis,"on"); end
    markerLineStyle = markerLineStyles(mod(index-1,numel(markerLineStyles))+1);
    plot(figAxis,xDataSets{index},yDataSets{index},markerLineStyle,...
        "DisplayName",legendNames(index));
end
hold(figAxis,"off");

if addLegend, legend(figAxis), end
end

%% --- Validation functions ---
function mustBeANamedParameterOrKeyFlag(paramName,resultsSets)
if paramName == QKDPlot.keyRateTag, return; end
if ~all(cellfun(@(x)ismember(paramName,string(fieldnames(x(1).currentParams))), resultsSets))
    throwAsCaller(MException("plotParameters:NotAParameter",...
        sprintf("%s was not a parameter name used in all the results sets or the key rate tag,'%s'.",paramName,QKDPlot.keyRateTag)))
end
end

function QKDSolverInputMustHaveExactlyOneScanParam(qkdInput)
numParams = cellfun(@numel,{qkdInput.scanParameters});
if any(numParams~=1)
    throw(MException("plotResults:MustHaveExactlyOneScanParam",...
        "The QKDSolverInput must have exactly one scan parameter for plotting."))
end
end
