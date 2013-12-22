function plotSOLProblemConvergence
% Plot convergence tests for SOL problem
    close all
    clear
    
    pathBase = '/Users/eshi/Research/gkeyllall/gkeyll-sims/eric/sol/solConvergence/';
    filenumEnd = 1;
    titles = {'Total Heat Flux', 'Ions + Electrons'};
    
    xSizes = {'16','32'};
    vSizes = {'32'};
    filenames1D = cell(1, length(xSizes)*length(vSizes));
    legendNames = cell(1, length(xSizes)*length(vSizes));
    dataName = 'heatFluxAtEdge';
 
    colors = {'r','b','m','k','c','g', 'y'};
    
    % Generate list of filenames
    for xIndex = 1:length(xSizes)
        for vIndex = 1:length(vSizes)
            legendNames{1, (xIndex-1)*length(vSizes) + vIndex} = [xSizes{xIndex},...
                'X',vSizes{vIndex},'V'];
            filenames1D{1, (xIndex-1)*length(vSizes) + vIndex} = [legendNames{1, (xIndex-1)*length(vSizes) + vIndex},...
                '_',dataName];
        end
    end
    
    
    % Plot all this data on one plot
    if ~isempty(filenames1D)
        % Cell to store all plot matrices
        plotCell1D = cell(1,size(filenames1D,2));
        for dataIndex = 1:size(filenames1D,2)
            % Get a matrix with x vals in col 1 and y vals in col 2
            plotCell1D{dataIndex} = get1DPlotMatrix(filenames1D{dataIndex}, filenumEnd, pathBase);
        end
        
        % Loop over total flux, ion flux, and electron flux
        for dataIndex = 2:3
            % Make first plot
            figure
            semilogx(plotCell1D{1}(:,1).*1e6,plotCell1D{1}(:,dataIndex),'r','LineWidth',1.5)
            % Loop over rest
            hold on
            if dataIndex == 2
                % Total only, every file
                for fileIndex = 2:size(filenames1D,2)
                    semilogx(plotCell1D{fileIndex}(:,1).*1e6,plotCell1D{fileIndex}(:,dataIndex),...
                        colors{mod(fileIndex-1,length(colors))+1},'LineWidth',1.5)
                end
            else
                % Ions + Electrons
                for fileIndex = 2:size(filenames1D,2)
                    semilogx(plotCell1D{fileIndex}(:,1).*1e6,plotCell1D{fileIndex}(:,dataIndex),...
                        colors{mod(fileIndex-1,length(colors))+1},'LineWidth',1.5)
                end
                
                semilogx(plotCell1D{1}(:,1).*1e6,plotCell1D{1}(:,dataIndex+1),'r','LineWidth',1.5)
                for fileIndex = 2:size(filenames1D,2)
                    semilogx(plotCell1D{fileIndex}(:,1).*1e6,plotCell1D{fileIndex}(:,dataIndex+1),...
                        colors{mod(fileIndex-1,length(colors))+1},'LineWidth',1.5)
                end
            end
            hold off
            
            ylim([0 5e9])
            xlim([0.01 1000])
            grid on
            legend(legendNames)
            xlabel('t (\mu s)')
            ylabel('Q (W/m^2)')
            title(titles{dataIndex-1})
        end
    end
end

function [functionVector, NDIM, basisDegree, vsNumCells] = getBasisInfo(dataname, pathBase)
    % For now just take first file, and num = 0
    filename = [pathBase,dataname,'_',num2str(0),'.h5'];
    infoStruct = h5info(filename,'/StructGridField');
    dataSize   = infoStruct.Dataspace.Size;
    NDIM       = length(h5readatt(filename,'/StructGrid','vsUpperBounds'));
    vsNumCells = double(h5readatt(filename,'/StructGrid','vsNumCells'));
    numBasis   = dataSize(1);
    % Figure out 2-D and 1-D basis functions
    if NDIM == 1
        basisDegree = numBasis - 1;
    elseif NDIM == 2
        switch numBasis
            case 4
                basisDegree = 1;
            case 8
                basisDegree = 2;
            case 12
                basisDegree = 3;
        end
    end

    functionVector = getSerendipityBasisFunctions(NDIM, basisDegree);
end

function plotMatrix = get1DPlotMatrix(dataname, filenumEnd, pathBase)
    plotCell = cell(filenumEnd+1,2);
    
    for filenum = 0:filenumEnd
        filename = [pathBase,dataname,'_',num2str(filenum),'.h5'];
%         infoStruct = h5info(filename,'/DataStruct/data');
%         infoStruct.Dataspace.Size
        timeVals = h5read(filename,'/DataStruct/timeMesh');
        dataVals = h5read(filename,'/DataStruct/data');
        plotCell{filenum+1,1} = timeVals';
        plotCell{filenum+1,2} = dataVals';
    end
    
    plotMatrix = cell2mat(plotCell);
end