function plotSOLProblemInitialConditions
% Plot one collision test for SOL problem
    close all
    clear
    
    pathBase = '/Users/eshi/Research/gkeyllall/gkeyll-tests/tests/slvrs/solInitialProfileTest/';
    filenumEnd = 10;
    titles = {'Total Heat Flux', 'Ions + Electrons'};
    
    filenames1D = {'uniform','nonuniform'};
    legendNames = {'uniform','non-uniform'};
    dataName = 'heatFluxAtEdge';
    
    filenames2D = {'numDensity','tempProfile'};
    subCellsPerLength = [4 4];
 
    colors = {'r','b','m','k','c','g', 'y'};
    
    % Plot 2-D Data first
        % Plot all this data on one plot
    if ~isempty(filenames2D)
        for dataSet = 1:size(filenames2D, 2)
            % Generate list of filenames
            filenamesToRead = cell(1,length(filenames1D));
            for fileIndex = 1:length(filenames1D)
                filenamesToRead{1, fileIndex} = [filenames1D{1, fileIndex},'_',filenames2D{1, dataSet}];
            end
            % Cell to store all plot matrices
            % Each row is all sets of data at a particular time
            plotCell2D = cell(filenumEnd+1, length(filenames1D));
            ndimList = zeros(1, length(filenames1D));
            vsNumCellsList = cell(1, length(filenames1D));
            for dataIndex = 1:size(filenamesToRead,2)
                % Get a matrix with x vals in col 1 and y vals in col 2
                [plotCell2D{dataIndex},ndimList(dataIndex),vsNumCellsList{dataIndex}] = get2DPlotMatrix(filenamesToRead{dataIndex}, filenumEnd, pathBase, subCellsPerLength);
            end
            
            % Plot both uniform and nonuniform data for this filename on
            % same plot
            figure
            hold on
            for dataIndex = 1:size(filenamesToRead,2)
                outputCell = plotCell2D{dataIndex};
                timeIndex = 0;
                
                outputMatrix = outputCell{timeIndex+1};
                NDIM = ndimList(dataIndex);
                vsNumCells = vsNumCellsList{dataIndex};

                if NDIM == 2
                    dataRowsPerX = vsNumCells(2)*subCellsPerLength(2);
                    vVals = outputMatrix( 1:dataRowsPerX,2);
                    fVals = outputMatrix( 1:dataRowsPerX,3);
                    plot(vVals, fVals, colors{mod(dataIndex-1,length(colors))+1},'LineWidth',1.5)
                elseif NDIM == 1
                    xVals = outputMatrix(:,1);
                    fVals = outputMatrix(:,2);
                    plot(xVals, fVals, colors{mod(dataIndex-1,length(colors))+1},'LineWidth',1.5)
                end            
            end
            hold off
            legend(legendNames,'Location','Best')
            title(filenames2D{dataSet})
        end
    end
    
    % Generate list of filenames
    for fileIndex = 1:length(filenames1D)
        filenames1D{1, fileIndex} = [filenames1D{1, fileIndex},'_',dataName];
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
                % Total only
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
                
                for fileIndex = 2:size(filenames1D,2)
                    semilogx(plotCell1D{fileIndex}(:,1).*1e6,plotCell1D{fileIndex}(:,dataIndex+1),...
                        colors{mod(fileIndex-1,length(colors))+1},'LineWidth',1.5)
                end
            end
            hold off
            
            ylim([0 5e9])
            xlim([0.01 1000])
            grid on
            legend(legendNames,'Location','Best')
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

function [plotCell, NDIM, vsNumCells] = get2DPlotMatrix(dataname, filenumEnd, pathBase, subCellsPerLength)
% Return a cell containing the full interpolated solution at a particular time
    [functionVector, NDIM, basisDegree, vsNumCells] = getBasisInfo(dataname, pathBase);
    
    plotCell = cell(1,filenumEnd+1);
    
    for filenum = 0:filenumEnd
        filename = [pathBase,dataname,'_',num2str(filenum),'.h5'];
        plotCell{filenum+1} = getInterpolatedPlotData(filename, NDIM, basisDegree, functionVector, subCellsPerLength);
    end
end