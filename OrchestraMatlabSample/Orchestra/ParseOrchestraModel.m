function [chemistryFilePath, varDefinitionTable, ioVariableList] = ParseOrchestraModel(modelDir)
    % Parse ORCHESTRA model files and prepare inputs for running model from
    % MATLAB

    %% Parse colum.dat file with basic varialbes and file names
    [fileNames, iniVarDefinitionTable] = ParseColumnDat([modelDir 'column.dat']);

    chemistryFilePath = ['/' modelDir fileNames.chemistry];

    %% Parse file with variables and their default values
    [varList, varVals] = ReadVariables([modelDir fileNames.input]);
    varDefinitionTable = [varList varVals];
    varDefinitionTable = vertcat(iniVarDefinitionTable, varDefinitionTable);
    
    %% Parse file with the list of IO variables
    [ioVariableList, ~] = ReadVariables([modelDir fileNames.outputformat]);
    
    %% End
    return

    
    
    %% Nested functions
    function [files, varsList] = ParseColumnDat(fileName)
        columnFile = textread(fileName, '%s', 'delimiter', '\n');    
        
        varsList = { };
        
        files = struct();
        
        for i = 1:numel(columnFile)
            line = columnFile{i};
            
            % Skip lines that are: empty, commments
            if (isempty(line) || (strcmp(line(1:2), '//')))
                continue
            end
            
            % Split strings to a separate words
            tk = TokenizeString(line);
            
            if (strcmp(tk{1}, '@Var:'))
                varsList = vertcat(varsList, [tk(2), str2double(tk(3))]);
            end
            
            if (strcmp(tk{1}, '@include:'))
                files.input = tk{2};
            end
            
            if (strcmp(tk{1}, '@include_global_vars:'))
                files.chemistry = [tk{2} '.inp'];
            end
            
            if (strcmp(tk{1}, '@include_vars:') && strcmp(tk{2}(end-3:end), '.txt'))
                files.outputformat = tk{2};
            end
        end
    end

    function [varList, varVals] = ReadVariables(fileName)
        varList = { };
        varVals = { };
        
        inpFile = textread(fileName, '%s', 'delimiter', '\n');
        
        for i = 1:numel(inpFile)
            line = inpFile{i};
            
            % Skip lines that are: empty, commments
            if (isempty(line) || (strcmp(line(1:2), '//')))
                continue
            end
            
            % Split strings to a separate words
            tk = TokenizeString(line);
            
            if (strcmp(tk{1}, 'Var:'))
                varList = tk(2:end)';
            end
            
            if (strcmp(tk{1}, 'Data:'))
                varVals = num2cell(str2double(tk(2:end)'));
            end
        end
    end

    function tokens = TokenizeString(inStr)
        [tokens, rem] = strtok(inStr);
        while ~isempty(rem)
            [tok, rem] = strtok(rem);
            if ~isempty(tok)
                tokens = [tokens, {tok}];
            end
        end
    end

end