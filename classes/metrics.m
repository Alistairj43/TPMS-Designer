classdef Metrics
    % Metrics        
    % The metrics class stores and evaluates various properties using other
    % classes within the TPMS Designer package
    %
    % Methods:
    %   M = metrics(); - create and returns an empty metrics object 'M'
    %
    %   M = fvMetrics(FV); - calculate and updates metrics based on the
    %         properties of a surfaceMesh FV
    %
    %   M = mechanicalMetrics(F);  - calcualte the mechanical metrics
    %           from a the solid-field object
    %
    %   out = export();  - export the metrics as a classless
    %         data structure
    %
    % Properties
    %     Computational:
    %         CPUtime, Nfaces, Nnodes, errorFlag
    %     Geometric:
    %         volume, volumeFraction, surfaceArea, relativeArea, rmsMC, rmsGC, thickness, poreDiameter  
    %     Z-Slice:
    %         thicknessAM, areaVar, areaMean, areaStd, areaMin, areaMax,
    %       LpbfSimple, LpbfRa_Max, LpbfError_Mean, LpbfError_Max, areaBelow30deg 
    %     Mechanical:
    %         elastic, poisson, shear, totalStiffness, zenerRatio
    %
    % TPMS Design Package - Metrics class
    % Created by Alistair Jones, RMIT University 2022.
    
    properties
        %Computational
        CPUtime
        Nfaces
        Nnodes
        errorFlag
        
        %Geometric
        meshVolume
        relativeVolume
        surfaceArea
        relativeArea
        rmsMC
        rmsGC
        thickness
        poreDiameter

        %Z-Slice Analysis Properties
        thicknessAM
        areaVar
        areaMean
        areaStd
        areaMin
        areaMax
        
        %Mechanical 
        elastic
        poisson
        shear
        totalStiffness
        zenerRatio
    end
    
    methods
        function M = Metrics()
            %Constructor
            if nargin==0
                return;
            end
        end
        
        function M = fvMetrics(M,FV)
            % Function to calculate metrics based on properties of a surfaceMesh.
            M.surfaceArea = FV.totalArea;
            M.meshVolume = FV.totalVolume;
            M.Nfaces = size(FV.faces,1);
            M.Nnodes = size(FV.vertices,1);
            try
                M.errorFlag = 0;
                if isfield(FV.Vproperty,'MC')
                    M.rmsMC = rms(FV.Vproperty.MC);
                    M.rmsGC = rms(FV.Vproperty.GC);
                end
            catch
                M.errorFlag = 'No Curvature';
            end
        end
        
        function M = mechanicalMetrics(M,F)
            % Function to calculate the basic mechanical metrics
            [sx, sy, sz] = size(F.property.solid);
            M.relativeVolume = sum(F.property.solid(1:sx-1,1:sy-1,1:sz-1),'all')/((sx-1)*(sy-1)*(sz-1));

            % Calculate Z-Slice metrics
            M.areaMean = mean(F.zSlices.area);
            M.areaStd = std(F.zSlices.area);
            M.areaVar = var(F.zSlices.area);
            M.areaMin = min(F.zSlices.area,[],'all');
            M.areaMax = max(F.zSlices.area,[],'all');
            M.thicknessAM = mean(F.zSlices.maxthickness);

            %Calculate pore diameter (Pad to account for periodicity)
            
            tempArray = bwdist(padarray(padarray(F.property.solid,ceil(F.res),'circular','both'),[1 1 1],1,'both')); 
            M.poreDiameter = 2*F.voxelSize(1)*max(tempArray,[],'all');

            %Calculate approximate wall thickness
            tempArray = bwdist(padarray(padarray(~F.property.solid,ceil(F.res),'circular','both'),[1 1 1],1,'both')); 
            M.thickness = 2*F.voxelSize(1)*max(tempArray,[],'all');

            %Calculate mechanical metrics if the stiffness tensor is full rank
            try
            if rank(F.CH)==6
                    S = inv(F.CH);
                    E = zeros(6,1);
                    E(1) = 1/S(1,1);
                    E(2) = 1/S(2,2);
                    E(3) = 1/S(3,3);
                    E(4) = 1/S(4,4);
                    E(5) = 1/S(5,5);
                    E(6) = 1/S(6,6);
                    M.elastic = E(1);
                    M.shear = mean(E(4:6));
                    M.poisson = -S(1,2)*E(1);
                    M.totalStiffness = M.elastic+2*M.shear*(1-M.poisson);
                    M.zenerRatio = 2*(1+M.poisson)*M.shear/M.elastic; %2*F.CH(4,4)/(F.CH(1,1)-F.CH(1,2));
            end
            catch
                    M.errorFlag = 'Invalid Tensor';
            end
        end
        
        function out = export(M)
            temp = string(fieldnames(M));
            for i = 1:length(temp)
                if isempty(M.(temp(i)))
                    out.(temp(i)) = NaN;
                else
                    out.(temp(i)) = M.(temp(i));
                end
            end
        end
        
    end
end

