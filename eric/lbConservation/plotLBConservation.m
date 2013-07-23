function plotLBConservation
% Plot data output from 1x1v LB operator
    close all
    clear
    
    subCellsPerLength = [4 4];
    pathBase = '/Users/dark1egion/Research/Gkeyll-Project/gkeyllall/gkeyll-tests/gkeyll-tests/tests/slvrs/lb-operator-test_';
    filenumEnd = 1;
    filenames2D = {'distf'};
    filenames1D = {'totalPtclEnergy'};
 
    colors = {'r','b','m','c','g','k', 'y'};
    legendNames = cell(1, filenumEnd+1);
    makeColor = false;
    
    % Plot all this data on one plot
    if ~isempty(filenames2D)
        % Cell to store all plot matrices
        % Each row is all sets of data at a particular time
        plotCell2D = cell(filenumEnd+1,size(filenames2D,2));
        ndimList = zeros(1,size(filenames2D,2));
        vsNumCellsList = cell(1,size(filenames2D,2));
        for dataIndex = 1:size(filenames2D,2)
            % Get a matrix with x vals in col 1 and y vals in col 2
            [plotCell2D{dataIndex},ndimList(dataIndex),vsNumCellsList{dataIndex}] = get2DPlotMatrix(filenames2D{dataIndex}, filenumEnd, pathBase, subCellsPerLength);
        end
        
        for dataIndex = 1:size(filenames2D,2)
            outputCell = plotCell2D{dataIndex};

            figure
            hold on
            for timeIndex = 0:filenumEnd
                outputMatrix = outputCell{timeIndex+1};
                NDIM = ndimList(dataIndex);
                vsNumCells = vsNumCellsList{dataIndex};

                if NDIM == 2
                    vVals = outputMatrix(1:vsNumCells(2)*subCellsPerLength(2),2);
                    fVals = outputMatrix(1:vsNumCells(2)*subCellsPerLength(2),3);
                    plot(vVals, fVals, colors{mod(timeIndex,length(colors))+1})
                elseif NDIM == 1
                    xVals = outputMatrix(:,1);
                    fVals = outputMatrix(:,2);
                    plot(xVals, fVals, colors{mod(timeIndex,length(colors))+1})
                end
                legendNames{timeIndex+1} = num2str(timeIndex);
            end

            hold off
            legend(legendNames)
            title(filenames2D{dataIndex})

            if NDIM == 2 && makeColor == true
                makeColorMap(outputMatrix, vsNumCells, subCellsPerLength)
            end
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
        
        figure
        hold on
        for dataIndex = 1:size(filenames1D,2)
            plot(plotCell1D{dataIndex}(:,1),plotCell1D{dataIndex}(:,2),colors{mod(dataIndex,length(colors))+1});
        end
        hold off
        legend(filenames1D)
        set(gca, 'YTickLabel', num2str(get(gca,'YTick')',10))
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
        tempTime = h5read(filename,'/DataStruct/timeMesh');
        tempData = h5read(filename,'/DataStruct/data');
        plotCell{filenum+1,1} = tempTime';
        plotCell{filenum+1,2} = tempData';
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

function makeColorMap(outputMatrix, vsNumCells, subCellsPerLength)
% NDIM = 2 assumed
    colorPlot = zeros(vsNumCells(2)*subCellsPerLength(2),vsNumCells(1)*subCellsPerLength(1));
    for rowIndex = 1:size(outputMatrix,1)
        colorPlot(outputMatrix(rowIndex,5),outputMatrix(rowIndex,4)) = outputMatrix(rowIndex,3);
    end
    
    figure
    imagesc(colorPlot)
    axis square
end