classdef bulkSize
    %Object to store and define gloal object data
    
    properties
        method
        bbox
        FV
    end
    
    methods
        function obj = bulkSize(method,data)
            %REGION Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                method string = "box";
                data = [];
            end
            obj.method = method;
            
            switch method
                case "box"
                    obj.bbox = [-data; data]./2;
                case "cylinder"
                    obj.bbox = [-data; data]./2;
                case "sphere"
                    obj.bbox = [-data -data -data; data data data];
                case "ellipsoid"
                    obj.bbox = [-data; data];
                case "FV"
                    obj.FV.faces = data.faces;
                    obj.FV.vertices = data.vertices;
                    obj.bbox = [min(obj.FV.vertices,[],1)-[1 1 1]; max(obj.FV.vertices,[],1)]+[1 1 1];
                otherwise % In case of invalid inputs
                    obj.bbox = [-1 -1 -1; 1 1 1]./2;
                    obj.method = "box";
            end
        end
    end
end

