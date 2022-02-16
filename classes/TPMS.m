classdef TPMS
    %TPMS is a class to work with TPMS-like structures
    %   TPMS is designed to work within the TPMS Designer toolbox
    %   TPMS can be parametrically defined using an equation, type,
    %   isovalue/volume fraction, resolution, cell size, number of cells,
    %   rotation system, and phase offset.
    %TPMS has methods to convert the implicit form into a surfaceMesh object, a
    %v3field object, and to populate a metrics object.
    %TPMS also accesses methods for visualising the resulting data
    %structures
    %
    % TPMS Design Package TPMS class
    % Created by Alistair Jones, RMIT University 2021.
    
    properties
        equation
        u
        type
        v
        vf
        res
        cellSize
        ncells
        Rxyz
        offset
        FV
        F
        M
    end
    
    methods
        function TPMS = TPMS(equation,type,v,vf,res,cellSize,ncells,Rxyz,offset)
            %CREATEDOE generates a DOE based so a single loop can be used for
            %generating a full-factorial experimental settup
            %   Detailed explanation goes here
            arguments
                equation = "Gyroid";
                type string = "Network";
                v = []; % isovalue
                vf = 0.5; % Volume fraction
                res (1,3) double = [30 30 30];
                cellSize (1,3) double = [10 10 10];
                ncells (1,3) double = [1 1 1];
                Rxyz (1,3) double = [0 0 0];
                offset (1,3) double = [0 0 0];
            end
            TPMS.equation = equation;
            TPMS.type = type;
            TPMS.vf = vf;
            TPMS.v = v;
            TPMS.res = res;
            TPMS.cellSize = cellSize;
            TPMS.ncells = ncells;
            TPMS.Rxyz =  Rxyz;
            TPMS.offset = offset;

            switch TPMS.equation
                case "Diamond"
                    TPMS.u = @(x,y,z) ((sin(x).*sin(y).*sin(z)+...
                        sin(x).*cos(y).*cos(z)+...
                        cos(x).*sin(y).*cos(z)+...
                        cos(x).*cos(y).*sin(z)));
                    A = 1/sqrt(2)*0.5844;
                case "Gyroid"
                    TPMS.u = @(x,y,z) (cos(x).*sin(y)+ ...
                        cos(y).*sin(z)+cos(z).*sin(x));
                    A = 2/3*0.4964;
                case "Primitive"
                    TPMS.u = @(x,y,z) (cos(x)+cos(y)+cos(z));
                    A = 1/3*0.8491;
                case "IWP"
                    TPMS.u = @(x,y,z) 2*(cos(x).*cos(y)+cos(y).*cos(z)+...
                        cos(z).*cos(x))...
                        -(cos(2*x)+ cos(2*y)+ cos(2*z));
                    A = 40.82;
                case "Neovius"
                    TPMS.u = @(x,y,z) (3*(cos(x)+cos(y)+cos(z))+4*cos(x).*cos(y).*cos(z));
                    A = 1/13;
                case "FRD"
                    TPMS.u = @(x,y,z) 4*cos(x).*cos(y).*cos(z)-(cos(2*x).*cos(2*y)...
                        +cos(2*y).*cos(2*z)+cos(2*z).*cos(2*x));
                    A = 40.82;
                case "Sinusoidal"
                    TPMS.u = @(x,y,z) sin(x)+...
                        sin(y)-(z-pi);
                    A = 1;
                case "Sphere"
                    TPMS.u = @(x,y,z)(x-pi).^2+(y-pi).^2+(z-pi).^2-2;
                    A = 1;
                case "P-normCube"
                    TPMS.u = @(x,y,z) (x-pi).^10+(y-pi).^10+(z-pi).^10;
                    A = 1;
                case "Taurus"
                    R = 1; r = 0.1; A = 1;
                    TPMS.u = @(x,y,z) (sqrt((x-pi).^2+(y-pi).^2)-R).^2+(z-pi).^2-r.^2;
                otherwise
                    TPMS.u = str2func(TPMS.equation);
                    A = 1;
            end

            if isempty(TPMS.v) % Convert from volume fraction to isovalue using linear scaling
                switch TPMS.type
                    case "network"
                        TPMS.v=(TPMS.vf-0.5)./A;
                    case "surface"
                        TPMS.v=(TPMS.vf)./(2*A);
                    case "single"
                        TPMS.v=TPMS.vf;
                end
            end
        end
        
        function out = export(TPMS)
            % Flattens structure to save only single value outputs
            out = TPMS.M.export;
            out.equation = TPMS.equation;
            out.type = TPMS.type;
            out.vf = TPMS.vf;
            out.v1 = max(TPMS.v,[],'all');
            out.v2 = min(TPMS.v,[],'all');
            out.res = TPMS.res;
            out.ax = TPMS.cellSize(1);
            out.ay = TPMS.cellSize(2);
            out.az = TPMS.cellSize(3);
            out.pitch = TPMS.Rxyz(1);
            out.roll = TPMS.Rxyz(2);
            out.yaw = TPMS.Rxyz(3);
            out.ox = TPMS.offset(1);
            out.oy = TPMS.offset(2);
            out.oz = TPMS.offset(3);
        end
        
        function h = plot(TPMS,ax,plottype,property1,property2,opts)
            % Handles the plotting, returning a handle to the figure
            arguments
                TPMS;
                ax = [];
                plottype string = "surfaceMesh";
                property1 string = "inclination Angle (deg)";
                property2 string = "gaussian curvature (mm^-^2)";
                opts = [];
            end
            
            if isempty(ax)
                figure; %Create a new figure and use it as the axis
                ax = gca;
            end
            
            switch plottype
                case "surfaceMesh"
                    h = TPMS.FV.plotMesh(property1,opts,ax);
                case "voxel"
                    h = TPMS.F.plotField(property1,"voxels",ax);
                case "orthoslice"
                    h = TPMS.F.plotField(property1,[],ax);
                case "slice"
                    n = 50;
                    h= TPMS.F.plotField(property1,n,ax);
                case "histogram"
                    h = plotHistogram(ax,TPMS.FV,property1,opts);
                case "histogram2"
                    h = plotHistogram2(ax,TPMS.FV,property1,property2,opts);
                case "polemap"
                    h = plotPoleFig(ax,TPMS.FV,3);
                case "tensor"
                    h = plotTensor(ax,TPMS.F.CH);
            end
        end
        
        
        function TPMS = properties(TPMS,usefield,mechanical,usetri,curvature)
            % Function to update all avaliable properties and metrics
            arguments
                TPMS;
                usefield = 1;
                mechanical = 1;
                usetri = 1;
                curvature = 'implicit';
            end
            tic;
            if isempty(TPMS.F)
                TPMS = TPMS.field;
            end
            if isempty(TPMS.M)
                TPMS.M = metrics;
            end
            if isempty(TPMS.FV)
                TPMS.FV = surfaceMesh('v3field',TPMS.F);
            end
            
            if usefield||mechanical
                TPMS.F = TPMS.F.calculateProperties; 
                if mechanical
                    TPMS.F = TPMS.F.homogenise;
                end
                TPMS.M = TPMS.M.mechanicalMetrics(TPMS.F);
            end

            if usetri||curvature
                TPMS.FV = TPMS.FV.calculateProperties(curvature,TPMS);
                TPMS.M = TPMS.M.fvMetrics(TPMS.FV);
                temp = TPMS.cellSize.*TPMS.ncells;
                TPMS.M.relativeArea = TPMS.M.surfaceArea./...
                    (2*(temp(1)*temp(2)+temp(2)*temp(3)+temp(3)*temp(1)));
            end
            
            TPMS.M.CPUtime = toc; %Calculate the total computational time
        end
        
        function TPMS = field(TPMS)  
            Ftemp = v3Field("empty",[],ceil(TPMS.res.*TPMS.ncells),0,0,[0 0 0],TPMS.cellSize.*TPMS.ncells);  % Move to discretized space

            % Resample if required
            vtemp = voxresize(TPMS.v,TPMS.res.*TPMS.ncells);
            
            %Apply transformations
            rZ = [cosd(TPMS.Rxyz(3)) -sind(TPMS.Rxyz(3)) 0; sind(TPMS.Rxyz(3)) cosd(TPMS.Rxyz(3)) 0; 0 0 1];
            rY = [cosd(TPMS.Rxyz(2)) 0 -sind(TPMS.Rxyz(2)); 0 1 0; sind(TPMS.Rxyz(2)) 0 cosd(TPMS.Rxyz(2))];
            rX = [1 0 0; 0 cosd(TPMS.Rxyz(1)) -sind(TPMS.Rxyz(1)); 0 sind(TPMS.Rxyz(1)) cosd(TPMS.Rxyz(1))];
            Axyz = [TPMS.cellSize(1)/(2*pi) 0 0; 0 TPMS.cellSize(2)/(2*pi) 0; 0 0 TPMS.cellSize(3)/(2*pi)];
            trans = [TPMS.offset(1) TPMS.offset(2) TPMS.offset(3)].*TPMS.cellSize;
            pad = [0; 0; 0]; T = Axyz*rZ*rY*rX;
            tform=affine3d([T pad; trans 1]);
            [X, Y, Z] = transformPointsInverse(tform,Ftemp.property.X,Ftemp.property.Y,Ftemp.property.Z);
            
            switch TPMS.type
                case "double"
                    vtemp2 = voxresize(TPMS.vf,TPMS.res.*TPMS.ncells);
                    Ftemp.property.U = (TPMS.u(X,Y,Z)-vtemp).*(TPMS.u(X,Y,Z)-vtemp2);
                case "surface"
                    Ftemp.property.U = (TPMS.u(X,Y,Z)-vtemp).*(TPMS.u(X,Y,Z)+vtemp);
                case "bounded"
                    Ftemp.property.U = (TPMS.u(X,Y,Z)-TPMS.v).*(TPMS.u(X,Y,Z)-TPMS.vf);
                otherwise
                    Ftemp.property.U = TPMS.u(X,Y,Z)-vtemp;
            end
            Ftemp.property.solid = double(Ftemp.property.U<0);
            TPMS.F = Ftemp;
        end
    end
end
