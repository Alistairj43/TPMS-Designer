classdef UnitCell
    %UnitCell is a class to work with cellular solid structures
    %   UnitCell is designed to work within the TPMS Designer toolbox A
    %   UnitCell can be parametrically defined using an equation, type, two
    %   geometry parameters, resolution, and cell size
    %
    % TPMS Design Package - UnitCell class Created by Alistair Jones, RMIT
    % University 2022.
    
    properties
        type % Type of Cell equation 
        equation % - Equation or path to file describing
        eqn % Equation showing parameters
        u % Data/Equation describing features
        v1 % Geometry Parameter 1
        v2 % Geometry Parameter 2
        res % Resolution (x,y,z)
        cellSize % Size of unit cell (x,y,z)
        FV % SurfaceMesh Object
        F % V3Field Object
        M % Metrics Object
        B % Object describing the bulk object
    end
    
    methods
        function UnitCell = UnitCell(type,equation,v1,v2,res,cellSize,B)
            %Constructor for to initialising the UnitCell object 
            %   Inputs:
            %       type - (string) The type of unit cell (default = "network")
            %       equation - (string) The equation defining the unit cell
            %           (default = "gyroid") 
            %       v1 - (double) First geometric parameter (defualt = 0.5) 
            %       v2 - (double) Second geometric parameter (defualt = 0.5) 
            %       res - (1,3)(double) Resolution of cell in x,y,z (default = [30 30 30]) 
            %       cellSize - (1,3)(double) Size of cell in x,y,z (default = [10 10 10])
            %       B - Structure descirbing the bulk size of the
            %           object, includes a minimum and maxium bounding box coord and possibly other fields to describe the region)
            %   Outputs:
            %   Unit Cell - (UnitCell) UnitCell object
            arguments
                type string = "network";
                equation = "Gyroid";
                v1 double = 0.5; % Geometry parameter 1
                v2 double = 0; % Geometry parameter 2
                res (1,3) double = [30 30 30];
                cellSize (1,3) double = [100 100 100];
                B = [];
            end
            UnitCell.equation = equation;
            UnitCell.type = type;
            UnitCell.v1 = v1;
            UnitCell.v2 = v2;
            UnitCell.res = res;
            UnitCell.cellSize = cellSize;
            UnitCell.M = metrics();

            if isempty(B) % Default bulk size is one unit cell starting at the origin
                UnitCell.B = bulkSize('box',cellSize);
            else
                UnitCell.B = B;
            end


            A = [];
            ax = 2*pi./cellSize(1); ay = 2*pi./cellSize(2); az = 2*pi./cellSize(3);
            switch UnitCell.equation
                case "Diamond"
                    UnitCell.u = @(x,y,z) ((sin(ax.*x).*sin(ay.*y).*sin(az.*z)+...
                        sin(ax.*x).*cos(ay.*y).*cos(az.*z)+...
                        cos(ax.*x).*sin(ay.*y).*cos(az.*z)+...
                        cos(ax.*x).*cos(ay.*y).*sin(az.*z)));
                    A = 1/sqrt(2)*0.5844;
                case "Gyroid"
                    UnitCell.u = @(x,y,z) (cos(ax.*x).*sin(ay.*y)+ ...
                        cos(ay.*y).*sin(az.*z)+cos(az.*z).*sin(ax.*x));
                    A = 2/3*0.4964;
                case "Primitive"
                    UnitCell.u = @(x,y,z) (cos(ax.*x)+cos(ay.*y)+cos(az.*z));
                    A = 1/3*0.8491;
                case "IWP"
                    UnitCell.u = @(x,y,z) 2*(cos(ax.*x).*cos(ay.*y)+cos(ay.*y).*cos(az.*z)+...
                        cos(az.*z).*cos(ax.*x))...
                        -(cos(2*x)+ cos(2*y)+ cos(2*z));
                    A = 40.82;
                case "Neovius"
                    UnitCell.u = @(x,y,z) (3*(cos(ax.*x)+cos(ay.*y)+cos(az.*z))+4*cos(ax.*x).*cos(ay.*y).*cos(az.*z));
                    A = 1/13;
                case "FRD"
                    UnitCell.u = @(x,y,z) 4*cos(ax.*x).*cos(ay.*y).*cos(az.*z)-(cos(2*x).*cos(2*y)...
                        +cos(2*y).*cos(2*z)+cos(2*z).*cos(2*x));
                    A = 40.82;
                case "Sinusoidal"
                    UnitCell.u = @(x,y,z) sin(ax.*x)+...
                        sin(ay.*y)-(z-pi);
                case "Sphere"
                    UnitCell.u = @(x,y,z)(x-pi).^2+(y-pi).^2+(z-pi).^2-2;
                case "P-normCube"
                    UnitCell.u = @(x,y,z) (x-pi).^10+(y-pi).^10+(z-pi).^10;
                case "Taurus"
                    R = 1; r = 0.1;
                    UnitCell.u = @(x,y,z) (sqrt((x-pi).^2+(y-pi).^2)-R).^2+(z-pi).^2-r.^2;
                case "BCC"
                    [UnitCell.u.nodes,UnitCell.u.struts] = readLattice('data/lattices/BCC.txt');
                case "BCCXYZ"
                    [UnitCell.u.nodes,UnitCell.u.struts] = readLattice('data/lattices/BCCXYZ.txt');
                otherwise
                    if endsWith(UnitCell.equation,'.txt')
                        [nodes,struts] = readLattice(UnitCell.equation);
                        UnitCell.u.nodes = nodes;
                        UnitCell.u.struts = struts;
                    else
                        UnitCell.u = str2func(UnitCell.equation);
                    end
            end

            
            if ~isempty(A)
                switch UnitCell.type
                    case "network"
                        UnitCell.v1 = -(UnitCell.v1-0.5)./A;
                    case "surface"
                        UnitCell.v1 = (UnitCell.v1)./(2*A);
                end
            end
        end
        
        function out = export(UnitCell)
            %Flatten structure to exportable values
            % Inputs:
            %   Unit Cell - (UnitCell) Self-referenced object
            % Outputs:
            %   out - Struct with flattened (numerical) outputs
            out = UnitCell.M.export;
            out.equation = UnitCell.equation;
            out.type = UnitCell.type;
            out.v1 = UnitCell.v1;
            out.v2 = UnitCell.v2;
            out.res = UnitCell.res(1);
            out.Lx = UnitCell.cellSize(1);
            out.Ly = UnitCell.cellSize(2);
            out.Lz = UnitCell.cellSize(3);
        end
        
        function h = plot(UnitCell,plottype,property1,property2,opts,ax)
            %Handles plotting returning a handle to the created figure
            % Inputs:
            %   Unit Cell - (UnitCell) Self-referenced object 
            %   plottype - (String) The type of plot Valid inputs include
            %       'SurfaceMesh' (defualt), 'voxel',
            %       'orthoslice', 'slice', 'histogram', 'histogram2',
            %       'polemap', 'tensor', 'lattice'
            %   property1 - (String) Name of the primary property being
            %       visualised. Defaults to 'inclination'
            %   property2 - (String) Name of the secondary property being
            %       visualised. Defaults to 'azimuth'. This is other used
            %       for vectormapping on a SurfaceMesh and histogram2 plots
            %   opts - Options, additional parameter based on plot type
            %       opts.fancy (boolean): use fancy graphics styling (default=1), 
            %       opts.caps (int): 0=dont plot, 1=plot without color, 2=plot with color  
            %       opts.s (int): select the slice height %
            %       (0-100) opts.n (int): set the number of bins for
            %           histogram/histogram2/polemap
            %   ax - graphics object to plot on
            % Outputs:
            %   h - Handle to the created graphical object
            arguments
                UnitCell;
                plottype string = "SurfaceMesh";
                property1 string = "inclination";
                property2 string = "azimuth";
                opts = [];
                ax = [];
            end
            

            switch plottype
                case "voxel"
                    h = UnitCell.F.plotField(property1,"voxel",opts,ax);
                case "orthoslice"
                    h = UnitCell.F.plotField(property1,"orthoslice",opts,ax);
                case "slice"
                    if isfield(opts,'slice')
                        s = opts.slice;
                    else
                        s = 50;
                    end
                    h= UnitCell.F.plotField(property1,s,ax);
                case "histogram"
                    h = plotHistogram(UnitCell,property1,opts,ax);
                case "histogram2"
                    h = plotHistogram2(UnitCell.FV,property1,property2,opts,ax);
                case "polemap"
                    if isfield(opts,'n')
                        n = opts.n;
                    else
                        n = 3;
                    end
                    h = plotPoleFig(UnitCell.FV,n,opts,ax);
                case "tensor"
                    h = plotTensor(UnitCell.F.CH,ax);
                case "lattice"
                    h = plotLattice(UnitCell.u.nodes,UnitCell.u.struts,ax);
                otherwise
                    h = UnitCell.FV.plotMesh(property1,property2,0.8,[],ax);
            end
        end
        
        
        function UnitCell = compute(UnitCell,computeCurvature,computeMechanical,computeMesh)
            %Function to generate the unit cell and compute various
            % properties Inputs:
            %   Unit Cell - (UnitCell) Self-referenced object
            %   computeCurvature - (String) The method for curvature
            %       including 'trimesh2', 'meyer2003', 'implicit', 'none'
            %       see SurfaceMesh.calculateProperties() for more info.
            %   computeMechanical - (logical) (default = 1) 
            %   computeMesh - (logical) (default = 1)
            % Outputs:
            %   Unit Cell - (UnitCell) Self-referenced object with
            %       calcualted properties
            arguments
                UnitCell;
                computeCurvature = 'trimesh2';
                computeMechanical = 1;
                computeMesh = 1;
            end

            % Start the timer
            tic; 


            % Compute SDF and properties
            if isempty(UnitCell.F)
                if strcmp(UnitCell.type,'lattice') %Handle Lattices
                    UnitCell.u.rstrut = UnitCell.v1;
                    UnitCell.u.rnode = UnitCell.v2;
                    UnitCell.F = v3Field('lattice',UnitCell.u,UnitCell.res,UnitCell.B);
                else
                    % Initialise the Field for the TPMS using equation u.
                    switch UnitCell.type
                        case "surface"
                            temp = @(x,y,z)(UnitCell.u(x,y,z)-UnitCell.v1+UnitCell.v2).*(UnitCell.u(x,y,z)+UnitCell.v1+UnitCell.v2);
                        case "double"
                            temp = @(x,y,z)(UnitCell.u(x,y,z)+UnitCell.v1).*(UnitCell.u(x,y,z)+UnitCell.v2);
                        case "single"
                            temp = @(x,y,z)UnitCell.u(x,y,z)+UnitCell.v1;
                        otherwise
                            temp = @(x,y,z)UnitCell.u(x,y,z)+UnitCell.v1;
                    end
                    UnitCell.F = v3Field("TPMS",temp,ceil(UnitCell.res),UnitCell.B);  % Move to discretized space centred at origin
                end
            end

            % Calculate derived properties
            UnitCell.F = UnitCell.F.calculateProperties(computeCurvature,UnitCell);

            % Generate the SurfaceMesh
            if computeMesh&&isempty(UnitCell.FV) 
                UnitCell.FV = surfaceMesh('V3Field',UnitCell.F);
                UnitCell.FV = UnitCell.FV.calculateProperties(computeCurvature,UnitCell);
                UnitCell.M = UnitCell.M.fvMetrics(UnitCell.FV);                
                temp = UnitCell.cellSize;
                UnitCell.M.relativeArea = UnitCell.M.surfaceArea./...
                    (2*(temp(1)*temp(2)+temp(2)*temp(3)+temp(3)*temp(1)));
            end

            % Run Homogenisation
            if computeMechanical  
                UnitCell.F = UnitCell.F.homogenise();
            end

            %Volume Metrics
            UnitCell.M = UnitCell.M.mechanicalMetrics(UnitCell.F); 

            %Calculate the total computational time
            UnitCell.M.CPUtime = toc; 
        end
    end
end
