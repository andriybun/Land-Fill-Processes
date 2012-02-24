function RunOrchestraCallSample()

    % First we have to add the directory containing OrchestraInterface.jar
    % file and other *.txt files with configuration to our MATLAB paths
    addpath('../Orchestra/');
    
    % Then specify the path to chemistry.inp
    chemistryFilePath = '/../Data/chemistry.inp';

    % We define the table of variables specifying a variable name and its
    % default value
    varDefinitionTable = {
        'Ca+2.liter', 1
        'H+.liter', -2.5
    	'CO2[g].tot', 11
        'I', 0.814
        };
    
    % Specify an array with the names of IO variables
    ioVariableList = {
        'H[Acetate].liter'
        'pH'
        };
    
    % Then specify the name or index of IO variable we want to have a loop over
    hAcetateIdx = 1;
    
    % And finally define the range of values for this variable
    hAcetateVal = 0.01:0.01:0.1;
    
    % Call Orchestra interface
    result = RunOrchestraInterface(chemistryFilePath, varDefinitionTable, ioVariableList, hAcetateIdx, hAcetateVal);
    
    disp(result);
end