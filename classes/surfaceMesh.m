classdef surfaceMesh
    %surfaceMesh is a class to work with triangulation mesh objects
    %   trimesh is designed to work within the TPMS Designer toolbox
    %   trimeshes may be initialised from an implicit/field or an stl file
    %   methods have been added for calculating properties and visualising
    %   the trimesh object
    %
    % TPMS Design Package trimesh class
    % Created by Alistair Jones, RMIT University 2021.
    
    properties
        name
        faces
        vertices
        facesC
        verticesC
        totalVolume
        totalArea
        Fproperty
        Vproperty
        color
    end
    
    methods
        function FV = surfaceMesh(varargin)
            %surfaceMesh Construct an instance of this class
            %   Creates a surfaceMesh object, and allows various
            %   callbacks to convert and calculate properties.
            % The first input is a string describing the data
            % FV = surfaceMesh() Creates an empty array
            % FV = surfaceMesh('stl','path/object.stl',..) imports .stl from path 
            % FV = surfaceMesh('fv',FV) directly specify FV structure with matching
            % FV = surfaceMesh('field',field) for use with v3Field
            %   properties
            % FV = surfaceMesh(...,0) calculate first derivative properties
            % FV = surfaceMesh(...,1) calculatre all properties (default)
            if nargin == 0
                return
            end
            
            % Convert input types
            [varargin{:}] = convertStringsToChars(varargin{:});            
            if contains(varargin{1},'fv','IgnoreCase',true)
                % Construct from existing FV structure
                fields = fieldnames(varargin{2});
                isfield(fields,properties(surfaceMesh));
                for i = 1:length(fields)
                    FV.(fields{i})= varargin{2}.(fields{i});
                end
                
                if ~isfield(fields,'name')
                    FV.name = inputname(2);
                end
                
            elseif contains(varargin{1},'stl','IgnoreCase',true)
                % Construct from stl import
                [FV.vertices, FV.faces, FV.color, ~] = readSTL(string(varargin(2)));
                FV.name = string(varargin(2));
                FV.vertices = FV.vertices-min(FV.vertices,[],1);
                
            elseif contains(varargin{1},'v3field','IgnoreCase',true)
                [temp, FV.vertices] = reducepatch(isosurface(varargin{2}.property.X,varargin{2}.property.Y,varargin{2}.property.Z,varargin{2}.property.U,0),1);
                if ~isempty(temp)
                    FV.faces(:,1) = temp(:,1); FV.faces(:,2) = temp(:,2); FV.faces(:,3) = temp(:,3);
                end
                [FV.facesC, FV.verticesC] = reducepatch(isocaps(varargin{2}.property.X,varargin{2}.property.Y,varargin{2}.property.Z,-varargin{2}.property.U,0,'above'),1); 
            end
            
            if ~isfield(FV,'color')
                FV.color = zeros(size(FV.faces,1),1);
            end
        end
        
        
        function [FV] = calculateProperties(FV,curvatureMethod,TPMS)
            arguments
                FV;
                curvatureMethod string = 'patch';
                TPMS = [];
            end
            %Calculate properties of a given structure
            %   Detailed explanation goes here
            %TRIANGULATIONPROPERTIES updates a triangulation(FV) object with more
            %properties per-face
            % Inputs:
            %   - data structure with fields FV.faces and FV.vertices
            % Outputs:
            %   - augmented data structure
            % Compute basic properties for cap regions (volume/area) 
            if ~isempty(FV.facesC)
                Fin = FV.facesC;
                Vin = FV.verticesC;
                d13= [(Vin(Fin(:,2),1)-Vin(Fin(:,3),1)), (Vin(Fin(:,2),2)-Vin(Fin(:,3),2)), (Vin(Fin(:,2),3)-Vin(Fin(:,3),3))];
                d12= [(Vin(Fin(:,1),1)-Vin(Fin(:,2),1)), (Vin(Fin(:,1),2)-Vin(Fin(:,2),2)), (Vin(Fin(:,1),3)-Vin(Fin(:,2),3))];
                cr = cross(d13,d12,2); % Edge 1-3x1-2 Cross-product
                crNorm = sqrt(cr(:,1).^2+cr(:,2).^2+cr(:,3).^2);
                AC = 0.5*sqrt(cr(:,1).^2+cr(:,2).^2+cr(:,3).^2);% Area of each triangle
                normalC = cr./crNorm;% Unit normal for each triangle
                centroidC = (Vin(Fin(:,1),:)+Vin(Fin(:,2),:)+Vin(Fin(:,3),:))/3;
                FV.totalArea = sum(AC,'all','omitnan');
                FV.totalVolume = sum(AC.*centroidC(:,3).*normalC(:,3),'all','omitnan'); % Volume flux/contribution of each triangle
            else
                FV.totalVolume = 0;
                FV.totalArea = 0;
            end
            
            % Compute properties for the main surface
            if ~isempty(FV.faces)
                vp = FV.Vproperty; p = FV.Fproperty;
                % Set up variables for use
                x = FV.vertices(:,1); y = FV.vertices(:,2);
                z = FV.vertices(:,3); tri = FV.faces;
                tri3d=triangulation(FV.faces,FV.vertices);
                ntemp = normalize(vertexNormal(tri3d),2,'norm');
                %vp.Nx = ntemp(:,1); vp.Ny = ntemp(:,2); vp.Nz = ntemp(:,3);
                vp.inclination = 180-real(acosd(ntemp(:,3)));                
                
                ntemp = normalize(faceNormal(tri3d),2,'norm');
                p.Nx = ntemp(:,1); p.Ny = ntemp(:,2); p.Nz = ntemp(:,3);
                [azi, inc, ~] = cart2sph(p.Nx,p.Ny,p.Nz);
                p.Finclination = 90+180/pi*inc;
                p.Fazimuth = 180/pi*azi;
                Fcentroid = incenter(tri3d);
                v1=[x(tri(:,2))-x(tri(:,1)),y(tri(:,2))-y(tri(:,1)),z(tri(:,2))-z(tri(:,1))];
                v2=[x(tri(:,3))-x(tri(:,1)),y(tri(:,3))-y(tri(:,1)),z(tri(:,3))-z(tri(:,1))];
                p.Farea = 0.5*vecnorm(cross(v1,v2,2),2,2);
                p.Fvolume = p.Farea.*Fcentroid(:,3).*p.Nz;
                FV.totalVolume = FV.totalVolume + sum(p.Fvolume,'all');
                FV.totalArea = FV.totalArea + sum(p.Farea,'all');
                
                % Basic Derived Properties
                vp.LPBFQuality = 0.5+0.5*tanh((vp.inclination-30)/20); 
                c1 = 2.204; c2 = 63.76; c3 = 0.06736; c4 = -7.843;
                vp.Ra = c1+c2./(1+exp(c3.*(vp.inclination+c4)));
                c1 = 11.8; c2 = 176.4; c3 = 0.08364; c4 = -25.87;
                vp.Rz = c1+c2./(1+exp(c3.*(vp.inclination+c4)));
                
                switch curvatureMethod
                    case 'trimesh'
                        % Trimesh curvature method (computationally
                        % efficient)
                        [vp.GC, vp.MC, vp.k1, vp.k2, vp.RMS, vp.ABS]= trimeshcurvature(x,y,z,tri3d);
                    case 'patch'
                        % Patch curvature method.
                        % more computationally expensive
                        useThird = 1;
                        [vp.GC, vp.MC, vp.k1, vp.k2, vp.RMS, vp.ABS] = patchcurvature(FV.faces,FV.vertices,useThird);
                    case 'implicit'
                        %Implicit curvature method - only valid for TPMS
                        if ~isdeployed
                            [vp.GC, vp.MC, vp.k1, vp.k2, vp.RMS, vp.ABS] = implicitcurvature(TPMS,FV.vertices);
                        end
                end
                
                if isfield(vp,'GC')
                    % 2nd Order Metrics - LBF manufacturability for mm
                    % curvature
                    c1 = 113.8; c2 = 260; c3 = 0.08119; c4 = 61.38; c5 = 180.6;
                    vp.LPBFerror = c1+c2./(1+exp(c3.*(vp.inclination-c4)+c5.*abs(vp.GC)));
                end
                
                FV.Fproperty = p;
                FV.Vproperty = vp;
            end
        end
 
        function h = plotMesh(FV,pName,ax,opt)
            % Function to plot a coloured mesh onto axis 'ax'
            % FV.plotMesh(ax)
            % FV.plotmesh(ax,property)
            %Defaults
            arguments
                FV;
                pName string = [];
                ax = [];
                opt = [];
            end
            h = []; % Create empty object to return
            if isempty(ax)
                %figure; ax = gca;
            end
            pID = convertStringsToChars(extractBefore(pName+" ", " "));
            if ~isempty(FV.faces)
                if isfield(FV,pID)
                    cData = FV.(pID);
                elseif isfield(FV.Fproperty,pID)
                    cData = FV.Fproperty.(pID);
                elseif isfield(FV.Vproperty,pID)
                    cData = FV.Vproperty.(pID);
                else
                    pName = [];
                    cData = [0.8 0.8 0.8];
                end
                if size(cData,1)==size(FV.vertices,1)
                    opts.FaceColor = 'interp';
                else
                    opts.FaceColor = 'flat';
                end
                
                if isfield(opt,'LineStyle')
                    opts.LineStyle = opt.LineStyle;
                end
                
                opts.LineWidth = 0.001;
                opts.EdgeAlpha = 0.1;
                
                h = patch(ax,'Faces',FV.faces,'Vertices',FV.vertices,'FaceVertexCData',cData,opts);
                xlabel(ax,"X"); ylabel(ax,"Y"); zlabel(ax,"Z");
                try 
                    c=colorbar(ax); c.Label.String = pName; colormap(ax,"jet");
                    ticks = unique(sort([caxisrange(1) get(c, 'YTick') caxisrange(2)]));
                    set(c, 'YTick', ticks);
                catch
                    % Error handling here
                end
            end
            
            if isfield(opt,'fancy')&&opt.fancy
                material(h,'shiny');
                view(ax,62,31); 
                lighting(ax,'gouraud');
                ax.XTick = []; ax.YTick = []; ax.ZTick = [];
                ax.BoxStyle = 'full';
                camlight(ax);
                axis(ax,'manual','vis3d','equal','tight');
            end
            
            if ~isempty(FV.facesC)
                alpha = 1;
                hold(ax,"on");
                h(2) = patch(ax,'Faces',FV.facesC,'Vertices',FV.verticesC,'FaceColor',[0.2,0.2,0.2],...
                    'FaceAlpha',alpha,'LineStyle','none');
                hold(ax,"off");
            end
            rotate3d(ax,'on');
        end

        
        function exportSTL(FV,filename)
            % Funtion to convert an FV structure to a Field
            % inputs:
            %   FV - surfaceMesh object
            %   filename - path/name to save file as
            % outputs: none
            arguments
                FV;
                filename string = "TPMSObject.stl";
            end
            if isempty(FV.facesC)
                tri=triangulation(FV.faces,FV.vertices);
            else
                tri=triangulation([FV.faces; FV.facesC+length(FV.vertices)],[FV.vertices; FV.verticesC]);
            end
            
            stlwrite(tri,filename,"binary");
        end

        function F = getField(FV, voxelSize)
            % Funtion to convert an FV structure to a Field
            % inputs:
            %   voxelSize - desired size of each voxel. default 1.0mm
            % outputs:
            %   F - v3Field object based on the distance-field of the surfaceMesh object
            arguments
                FV;
                voxelSize (1,1) double = 1.0;
            end
            [lower, upper] = bounds(FV.vertices);
            res = ceil((upper-lower)./voxelSize);
            F = v3Field("FV",FV,res,0,0,lower,upper);
            F.name = FV.name;
        end
    end
end

