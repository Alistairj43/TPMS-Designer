classdef metrics
    %metrics        
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
    %           from a field-based representation "F" - (requires homogenisation).
    %
    %   out = export();  - export the metrics as a classless
    %         data structure
    %
    % Properties
    %     Mechanical Properties/Moduli:
    %         elastic, poisson, shear, zenerRatio, bulk, total
    %     Physical Properties:
    %         volume, volumeFraction, surfaceArea, rmsMC, areaBelow30deg 
    %     L-LBF Properties
    %         LpbfSimple, LpbfRa_Max, LpbfError_Mean, LpbfError_Max
    %     General Properties
    %         CPUtime, errorFlag
    %
    % TPMS Design Package metrics class
    % Created by Alistair Jones, RMIT University 2021.
    
    properties
        CPUtime
        Nfaces
        Nnodes
        
        volume
        volumeFraction
        surfaceArea
        relativeArea
        rmsMC
        rmsGC
        thickness
        poreDiameter
        areaVar
        areaMean
        areaStd
        areaMin
        areaMax
        
        elastic
        poisson
        shear
        totalStiffness
        zenerRatio
        errorFlag
    end
    
    methods
        function M = metrics()
            %Constructor
            if nargin==0
                return;
            end
        end
        
        function M = fvMetrics(M,FV)
            % Function to calculate metrics based on properties of a surfaceMesh.
            M.surfaceArea = FV.totalArea;
            M.volume = abs(FV.totalVolume);
            M.Nfaces = size(FV.faces,1);
            M.Nnodes = size(FV.vertices,1);
            try
                M.errorFlag = 0;
                if isfield(FV.Vproperty,'MC')
                    M.rmsMC = rms(FV.Vproperty.MC);
                    M.rmsGC = rms(FV.Vproperty.GC);
                end
            catch
                M.errorFlag = 1;
            end
        end
        
        function M = mechanicalMetrics(M,F)
            % Function to calculate the basic mechanical metrics
            [sx, sy, sz] = size(F.property.solid);
            M.volumeFraction = sum(F.property.solid(1:sx-1,1:sy-1,1:sz-1),'all')/((sx-1)*(sy-1)*(sz-1));
            if ~isempty(F.CH)
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
                M.zenerRatio = 2*(1+M.poisson)*M.shear/M.elastic;%2*F.CH(4,4)/(F.CH(1,1)-F.CH(1,2));
            end
            M.areaMean = mean(F.zSlices.area);
            M.areaStd = std(F.zSlices.area);
            M.areaVar = var(F.zSlices.area);
            M.areaMin = min(F.zSlices.area,[],'all');
            M.areaMax = max(F.zSlices.area,[],'all');
            M.thickness = mean(F.zSlices.thickness);
            temp = bwdist(padarray(F.property.solid,F.res,'circular'));
            M.poreDiameter = (F.xq(2)-F.xq(1))*max(temp,[],'all');
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

