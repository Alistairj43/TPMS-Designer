classdef v3Field
    %v3Field is a class to deal with 3D volume data
    %   v3Fields can be created from an implicit function, trimesh, or
    %   lattice structure (defined by the nodes and struts).
    %
    %   The v3Field class comes with visualisation options to look at a
    %   particular slice, orthographic slice viewer, or voxelated object
    %
    %   v3Fields representing periodic unit cells may be assessed 
    %   for mechanical properties using homogenisation with periodic boundary conditions 
    %
    % TPMS Design Package v3Field class
    % Created by Alistair Jones, RMIT University 2021.
    
    properties
        name
        lower
        upper
        res
        xq
        yq
        zq
        property
        zSlices
        CH
    end
    
    methods
        function F = v3Field(method,data,res,rstrut,rnode,lower,upper)
            %FIELD Construct an instance of this class
            %     optional inputs:
            %             res (3,1) double = [30; 30; 30];
            %             lower (3,1) double = [0; 0; 0];
            %             upper (3,1) double = [10; 10; 10];
            %             method -> defines the method to be used = "flow"
            %             data -> the data to generate the volume field
            %             from
            arguments
                method string = "sample";
                data = [];
                res (1,3) double = [30 30 30];
                rstrut (1,1) double = 0.05;
                rnode (1,1) double = 0.05;
                lower (1,3) double = [0 0 0];
                upper (1,3) double = [1 1 1];
            end
            F.res = res;
            F.lower = lower;
            F.upper = upper;
            F.xq = linspace(lower(1),upper(1),res(1));
            F.yq = linspace(lower(2),upper(2),res(2));
            F.zq = linspace(lower(3),upper(3),res(3));
            [F.property.X,F.property.Y,F.property.Z] = meshgrid(F.xq,F.yq,F.zq);
            switch method
                case "FV"
                    temp.faces = data.faces;
                    temp.vertices = 1+(res-1).*((data.vertices-lower)./(upper-lower));
                    F.property.surface = 1.0*polygon2voxel(temp,[res(2),res(1),res(3)],'wrap',true);
                    F.property.solid = imfill(1.0*F.property.surface);
                    F.property.U = double(1.0.*(bwdist(F.property.solid)-bwdist(1.0-F.property.solid)));
                case "data"
                    F.property.U = data;
                    F.property.solid = 1.0*(F.property.U<=0);
                case "lattice"
                    if isstring(data)
                        filepath = strcat("data/lattices/",data,".txt");
                    else
                        filepath = "data/lattices/BCC.txt";
                    end
                    [F.property.solid] = voxelateLattice(F.property.X,F.property.Y,F.property.Z,[F.lower; F.upper],filepath,rstrut,rnode);
                    F.property.U =  double(1.0.*(bwdist(F.property.solid)-bwdist(1.0-F.property.solid)));
            end
        end
        
        function F = calculateProperties(F,method,TPMS)
            %Evaluate the values of the field
            %Inputs
            %   method - method for volumetric differentiation
            %   TPMS - the TPMS object for implicit differentiation (optional)
            arguments
                F;
                method string = "sobel";
                TPMS = [];
            end
            p = F.property;

            if strcmp("implicit",method)&&~isempty(TPMS) % Use implicit differentiation to calculate proeprties
                
            else % Use field integration to calculate properties
                [~,p.azimuth,elevation] = imgradient3(imgaussfilt3(F.property.U, 0.5),method);
                p.inclination = 90-elevation;
            end

            
            % Manufacturability using elipsoidal weighting filter -> 1 = perfect, 0 = empty space
            n = 9; nq = linspace(-1,1,n); [q1, q2, q3] = ndgrid(nq,nq,nq);
            c = 2.5; c1 = 0; c2 = 0.5;
            H = 1/(c^3*(2*pi)^(3/2)).*exp(-(q1.^2+q2.^2+q3.^2)/(2*c^2)); % Make weightings with gaussian distribution
            H(:,:,ceil((n/2)+1):n) = c1.*H(:,:,ceil((n/2)+1):n); % Set top half equal to 0.1
            H(:,:,ceil((n/2))) = c2.*H(:,:,ceil((n/2))); % Set current plane properties
            
            
            temp = zeros(size(F.property.solid,1)+2*n,size(F.property.solid,2)+2*n,size(F.property.solid,3)+n); % "Air" Padding
            temp(:,:,1:n) = ones(size(F.property.solid,1)+2*n,size(F.property.solid,2)+2*n,n); % "Build Platten" Padding Below
            temp(n+1:n+size(F.property.solid,1),n+1:n+size(F.property.solid,2),n+1:n+size(F.property.solid,3)) = F.property.solid; % Placing the Part
            temp = min(1,imfilter(temp,H)); % Apply convolutional filter
            p.manufacturability = F.property.solid-1+temp(n+1:n+size(F.property.solid,1),n+1:n+size(F.property.solid,2),n+1:n+size(F.property.solid,3)).*F.property.solid;
            F.property = p;
            for i = 1:size(F.property.solid,3)
                % Pad to account for periodic boundary conditions
                temp = padarray(F.property.solid,[1 1],0);
                XY.thickness(i) = max(bwdist(~temp),[],'all').*mean((F.upper-F.lower)./F.res);
                XY.area(i) = sum(F.property.solid(:,:,i),'all')./(F.res(1)*F.res(2));
            end
            F.zSlices = XY;
        end
        
        function F = homogenise(F,E1,v1,E2,v2)
            %homogenise the unit cell using periodic boundary conditions
            arguments
                F;
                E1 (1,1) double = 1; % Elastic Modulus solid
                v1 (1,1) double = 0.33; % Poisson's Ratio solid
                E2 (1,1) double = 0; % Elastic Modulus secondary
                v2 (1,1) double = 0; % Poisson's Ratio secondary
            end

            lambda1 = E1*v1/((1+v1)*(1-2*v1));
            mu1 = E1/(2+2*v1);
            lambda2 = E2*v2/((1+v2)*(1-2*v2));
            mu2 = E2/(2+2*v2);
            try
                F.CH = homo3D(F.upper(1),F.upper(2),F.upper(3),lambda1,mu1,lambda2,mu2,F.solid);
            catch % If homogenisation fails return NaNs
                F.CH = NaN(6,6);
            end
        end
        
        
        function h = plotField(F,pName,zslice,ax)
            %% Uses orthosliceviewer or 'slice' to visualise a v3field
            % Inputs:
            %   F - (self)
            %   ax - Axis to plot field on
            %   pName - property to investigate
            %   zslice - z height as a percentage of the total height to
            %           slice (between 0 and 100)
            % Outputs:
            %   h - handle to created object
            % Example:
            %   h = F.plotField('Nz',50)
            %   Credit: Alistair Jones, 12/11/2020, RMIT University
            arguments
                F;
                pName string = 'U';
                zslice = [];
                ax = [];
            end
            
            if isempty(ax) % If not axis is provided, create one
                figure;
                ax = gca;
            end
            
            pID = convertStringsToChars(extractBefore(pName+" ", " "));
            if isfield(F.property,pID)
                cData = F.property.(pID);
            else
                cData = F.(pID);
            end
                         
            if isempty(zslice)
                %Orthoslice
                h = orthosliceViewer(cData,'Parent',ax.Parent,'Colormap',[0 0 0; flipud(jet)]);
            elseif isnumeric(zslice)
                %Slice on axis
                z = F.lower(3)+zslice/100*range(F.zq,'all');
                [Xq, Yq, Zq] = ndgrid(F.xq,F.yq,z);
                xyslice = interp3(F.Y,F.X,F.Z,cData,Xq,Yq,Zq);
                h = surf(ax,Yq,Xq,Zq,xyslice,'LineStyle','none');
                c = colorbar; c.Label.String = pName; colormap(ax,flipud(jet));
                caxis(ax,'manual');
            else
                %Voxelated
                cData = cData+F.solid*1i;
                [h,~,~,~,~,~,~] = voxelSurf(cData,false);
            end
        end
    end
end

