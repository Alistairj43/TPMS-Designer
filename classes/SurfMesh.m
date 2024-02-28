classdef SurfMesh
    %SurfMesh is a class to work with triangulation mesh objects
    %   designed to work within the TPMS Designer toolbox trimeshes may be 
    %   initialised from an implicit/field or an stl file 
    %   methods have been added for calculating properties and visualising
    %
    % TPMS Design Package - SurfMesh class
    % Created by Alistair Jones, RMIT University 2022.
    
    properties
        name
        faces
        vertices
        facesc
        verticesc
        totalVolume
        totalArea
        Vproperty % Vertex property
        Fproperty % Face Property
        Qproperty % Vector Properties
    end
    
    methods
        function FV = SurfMesh(varargin)
            % SurfMesh Construct an instance of this class
            %   Creates a SurfMesh object, and allows various
            %   callbacks to convert and calculate properties.
            % Inputs: varargin
            %   FV = SurfMesh() Creates an empty array
            %   FV = SurfMesh('stl','path/object.stl',..) imports .stl from path 
            %   FV = SurfMesh('fv',FV) directly specify FV structure with matching
            %   FV = SurfMesh('TPMS',UnitCell) for use with a TPMS
            %   FV = SurfMesh(...,0) calculate first derivative properties
            %   FV = SurfMesh(...,1) calculatre all properties (default)
            if nargin == 0
                return
            end
            
            % Convert input types
            [varargin{:}] = convertStringsToChars(varargin{:});            
            if contains(varargin{1},'fv','IgnoreCase',true)
                % Construct from existing FV structure
                fields = fieldnames(varargin{2});
                isfield(fields,properties(SurfMesh));
                for i = 1:length(fields)
                    FV.(fields{i})= varargin{2}.(fields{i});
                end
                
                if ~isfield(fields,'name')
                    FV.name = inputname(2);
                end
            elseif contains(varargin{1},'stl','IgnoreCase',true)
                % Construct from stl import
                TR = stlread(string(varargin(2)));
                FV.faces = TR.ConnectivityList;
                FV.vertices = TR.Points;
                FV.name = string(varargin(2)); 
            elseif contains(varargin{1},'TPMS','IgnoreCase',true)
                [f, v] = isosurface(varargin{2}.F.property.X,varargin{2}.F.property.Y,varargin{2}.F.property.Z,varargin{2}.F.property.U,0);
                [FV.faces,FV.vertices]=cleanMesh(f, v);
                [fc, vc] = isocaps(varargin{2}.F.property.X,varargin{2}.F.property.Y,varargin{2}.F.property.Z,-varargin{2}.F.property.U,0,'above');
                if ~isempty(fc)
                    fc = fliplr(fc); % Fix normals to point outwards
                    [FV.facesc, FV.verticesc] = cleanMesh(fc, vc);
                end
            end
        end
        
        
        function [FV] = calculateProperties(FV,curvatureMethod,TPMS)
            arguments
                FV;
                curvatureMethod string = 'patch';
                TPMS = [];
            end
            %Calculate properties of a given structure
            %   Inlcuding curvatures, surface normals and derived
            %   proeprties
            % Inputs:
            %   FV - SurfMesh object
            %   curvatureMethod - (string) Describing method that should 
            %       be used to calculate curvature (default = 'trimesh2')
            %   TPMS - UnitCell object used for implicit differentiation
            %       (only for TPMS-type structures)
            % Outputs:
            %   FV - SurfMesh object with additional properties
            
            F = FV.faces;
            V = FV.vertices;
            
            % Compute properties for the main surface
            if size(V,2)==3&&size(F,2)==3
                p = FV.Vproperty; fp = FV.Fproperty;
                
                % Set up variables for use
                x = V(:,1); y = V(:,2); z = V(:,3);

                % Calculate normals
                FaceNormals=cross(V(F(:,3),:)-V(F(:,2),:),V(F(:,1),:)-V(F(:,3),:));
                FaceNormals=normr(FaceNormals); % Normalise to unit normals
                fp.Zheight = mean(z(F),2);
                fp.Nx = FaceNormals(:,1); fp.Ny = FaceNormals(:,2); fp.Nz = FaceNormals(:,3);
                [qp.UnitNormal,Avertex,Acorner,up,vp]=computeVertexNormals(FV,FaceNormals);
                p.inclination = 180-real(acosd(qp.UnitNormal(:,3))); 
                p.Nx = qp.UnitNormal(:,1); p.Ny = qp.UnitNormal(:,2); p.Nz = qp.UnitNormal(:,3);

                % Global Mesh Properties
                v1=[x(F(:,2))-x(F(:,1)),y(F(:,2))-y(F(:,1)),z(F(:,2))-z(F(:,1))];
                v2=[x(F(:,3))-x(F(:,1)),y(F(:,3))-y(F(:,1)),z(F(:,3))-z(F(:,1))];
                fp.Farea = 0.5*vecnorm(cross(v1,v2,2),2,2);
                fp.Fvolume = fp.Farea.*fp.Zheight.*fp.Nz;
                FV.totalVolume = sum(fp.Fvolume,'all','omitnan');
                FV.totalArea = sum(fp.Farea,'all','omitnan');

                if ~isempty(FV.facesc)

                end

                % Spherical Coordiante Transform
                [azi, inc, ~] = cart2sph(fp.Nx,fp.Ny,fp.Nz);
                fp.Finclination = 90+180/pi*inc;
                fp.Fazimuth = 180/pi*azi;
                [azi, ~, ~] = cart2sph(p.Nx,p.Ny,p.Nz);
                p.azimuth = 180/pi*azi;
                
                % Curvature Properties
                switch curvatureMethod
                    case {'trimesh2','Trimesh2'}
                        % Based on: "Estimating Curvatures and Their Derivatives on Triangle Meshes" 
                        % by Szymon Rusinkiewicz (2004) and according to its C implementation trimesh2 (More accurate)
                        [p.k1,p.k2,kDir1,kDir2,p.wfp]=curvatureTrimesh2(FV,qp.UnitNormal,FaceNormals,Avertex,Acorner,up,vp);
                        qp.kDir1 = kDir1.*p.k1;
                        qp.kDir2 = kDir2.*p.k2;
                    case {'meyer2003','Meyer2003'}
                        % Based on: "Discrete differential-geometry operators for triangulated 2-manifolds." 
                        % by Meyer, M., Desbrun, M., Schröder, P., & Barr, A. H. (2003).
                        [p.GC, p.MC, p.k1, p.k2]= curvatureMeyer2003(x,y,z,F);
                    case {'implicit','Implicit'}
                        % Implicit differentiation curvature method
                        % (accurate and fast when the equation is
                        % known/simple
                        [p.GC, p.MC, p.k1, p.k2] = curvatureImplicit(TPMS.u,FV.vertices,TPMS.tform);
                    case {'kroon','patch','Kroon'}
                        % Use Kroon Method using 3rd order vertices
                        [p.GC, p.MC, p.k1, p.k2, ~, ~]=curvatureKroon(F,V,1);
                    otherwise
                        p.k1 = zeros(size(FV.vertices,1),1);
                        p.k2 = p.k1;
                end
                p.GC = p.k1.*p.k2;
                p.MC = (p.k2+p.k1)/2;
                p.ABS = abs(p.k1)+abs(p.k2);
                p.RMS = real(sqrt(4*p.MC.^2-2*p.GC));


                % Basic Derived Properties
                p.LPBFQuality = 0.5+0.5*tanh((p.inclination-30)/20); 
                c1 = 2.204; c2 = 63.76; c3 = 0.06736; c4 = -7.843;
                p.Ra = c1+c2./(1+exp(c3.*(p.inclination+c4)));
                c1 = 11.8; c2 = 176.4; c3 = 0.08364; c4 = -25.87;
                p.Rz = c1+c2./(1+exp(c3.*(p.inclination+c4)));

                % Add Any Other Per-Face/Per Vertex Calculated Values Here
                if isfield(p,'GC')
                    % 2nd Order Metrics - LBF manufacturability for mm
                    % curvature
                    c1 = 113.8; c2 = 260; c3 = 0.08119; c4 = 61.38; c5 = 180.6;
                    p.LPBFerror = c1+c2./(1+exp(c3.*(p.inclination-c4)+c5.*abs(p.GC)));
                end

                FV.Fproperty = fp;
                FV.Vproperty = p;
                FV.Qproperty = qp;
            end
        end

        function h = plotMesh(FV,pName,p2Name,caps,opts,ax)
            % Function to plot a coloured mesh onto axis 'ax'
            %Defaults
            arguments
                FV SurfMesh;
                pName string = [];
                p2Name string = [];
                caps double = 0.6;
                opts = [];
                ax = [];
            end
            
            h = []; % Create empty object to return
            if isempty(ax)
                figure; ax = gca; rotate3d(ax,'on'); view(ax,62,31);
            end

            pID = convertStringsToChars(extractBefore(pName+" ", " "));
            if ~isempty(FV.vertices)
                if isfield(FV,pID)
                    cData = FV.(pID);
                elseif isfield(FV.Fproperty,pID)
                    cData = FV.Fproperty.(pID);
                elseif isfield(FV.Vproperty,pID)
                    cData = FV.Vproperty.(pID);
                else
                    pName = 'Surface Region';
                    cData = zeros(size(FV.faces,1),1);
                end
                

                opts.AlphaDataMapping = 'none'; opts.FaceColor = 'flat';
                opts.LineWidth = 0.001; opts.EdgeAlpha = 0.1;
                
                h = patch(ax,'Faces',FV.faces,'Vertices',FV.vertices,'FaceVertexCData',cData,opts);
                xlabel(ax,"X"); ylabel(ax,"Y"); zlabel(ax,"Z");
                c=colorbar(ax); c.Label.String = pName; colormap(ax,"jet");
                lighting(ax,'flat');
                ax.BoxStyle = 'full';
                axis(ax,'manual','vis3d','equal','tight');

                %Add function to plot vectors kDir1, kDir2, VertexNormals, FaceNormals%
                pID2 = convertStringsToChars(extractBefore(p2Name+" ", " "));
                if isfield(FV.Qproperty,pID2)
                    qData = FV.Qproperty.(pID2);
                    hold(ax,"on");
                    quiver3(ax,FV.vertices(:,1),FV.vertices(:,2),FV.vertices(:,3),qData(:,1),qData(:,2),qData(:,3),2);
                    hold(ax,"off");
                end
            end


            if caps&&~isempty(FV.facesc)
                opts.FaceAlpha = caps; opts.LineStyle = 'none'; hold(ax,'on');
                h(2) = patch(ax,'Faces',FV.facesc,'Vertices',FV.verticesc,opts,'FaceColor','black'); hold(ax,'off');
            end
        end
        
        function exportSTL(FV,filename)
            % Funtion to convert an FV structure to a Field
            % inputs:
            %   FV - SurfMesh object
            %   filename - path/name to save file as
            % outputs: none
            arguments
                FV;
                filename string = "TPMSObject.stl";
            end
            
            % Merge isocap with isosurface
            try
                F = [FV.faces; FV.facesc+length(FV.vertices)];
                V = [FV.vertices; FV.verticesc];
            catch
                F = FV.faces;
                V = FV.vertices;
            end

            tri=triangulation(F,V); % Convert to triangulation object
            stlwrite(tri,filename,"binary"); % Write using matlab default function
        end

        function F = getField(FV, voxelSize)
            % Funtion to convert an FV structure to a Field
            % inputs:
            %   voxelSize - desired size of each voxel. default 1.0mm
            % outputs:
            %   F - v3Field object based on the distance-field of the SurfMesh object
            arguments
                FV;
                voxelSize (1,1) double = 1.0;
            end
            region = bulkSize('FV',FV);
            F = v3Field("FV",FV,voxelSize,region); %F = v3Field(method,data,res,region,tform)
            F.name = FV.name;
        end

        function FEModelOut = computeFEMesh(FV, meshSize, meshGrowthRate, meshOrder)
            % Funtion to convert an FV structure to a Field
            % inputs:
            %   meshSize - desired size of mesh elements
            %   meshGrowthRate - grwoth rate
            %   meshOrder - order of elements ('linear or quadratic')
            % outputs:
            %   FEMESH - FEMesh object from the pde toolbox
            arguments
                FV;
                meshSize = 0.1;
                meshGrowthRate=1.1;
                meshOrder='quadratic';
            end

            FV.exportSTL("tempMesh.stl"); % Export temporarily
            FEModelOut = createpde(1); % Create model
            importGeometry(FEModelOut,"tempMesh.stl"); % Assign geometry
            generateMesh(FEModelOut); % Generate the mesh
        end
    end
end

