classdef v3Field
    %v3Field is a class to deal with 3D volume data
    %   v3Fields can be created from an implicit function, surfaceMesh, or
    %   lattice structure (defined by the nodes and struts).
    %
    %   The v3Field class comes with visualisation options to look at a
    %   particular slice, orthographic slice viewer, or voxelated object
    %
    %   v3Fields representing periodic unit cells may be assessed 
    %   for mechanical properties using homogenisation with periodic boundary conditions 
    %
    % TPMS Design Package v3Field class
    % Created by Alistair Jones, RMIT University 2022.
    
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
            %v3Field(method,data,res,rstrut,rnode,lower,upper) construct an instance of this class
            %   optional inputs:
            %             res (3,1) double = [30; 30; 30];
            %             lower (3,1) double = [0; 0; 0];
            %             upper (3,1) double = [10; 10; 10];
            %             method -> defines the method to be used = "flow"
            %             data -> the data to generate the v3field
            %   returns: the v3Field object
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
        
        function F = homogenise(F,E1,v1,E2,v2,method)
            %homogenise the unit cell using periodic boundary conditions
            arguments
                F;
                E1 (1,1) double = 1; % Elastic Modulus solid
                v1 (1,1) double = 0.33; % Poisson's Ratio solid
                E2 (1,1) double = 0; % Elastic Modulus secondary
                v2 (1,1) double = 0; % Poisson's Ratio secondary
                method string = 'pcg'; % Solver method
            end

            lambda1 = E1*v1/((1+v1)*(1-2*v1));
            mu1 = E1/(2+2*v1);
            lambda2 = E2*v2/((1+v2)*(1-2*v2));
            mu2 = E2/(2+2*v2);
            try
                F.CH = homo3D(F.upper(1),F.upper(2),F.upper(3),lambda1,mu1,lambda2,mu2,...
                    F.property.solid(1:(F.res(1)-1),1:(F.res(2)-1),1:(F.res(3)-1)),method);
            catch % If homogenisation fails return NaN
                F.CH = NaN(6,6);
            end
        end
        
        
        function h = plotField(F,pName,zslice,ax,opt)
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
                ax = gca;
                opt = [];
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
                xyslice = interp3(F.xq,F.yq,F.zq,cData,Xq,Yq,Zq);
                h = surf(ax,Yq,Xq,Zq,xyslice,'LineStyle','none');
                colormap(ax,"jet");
            else
                h = plotVoxel(ax,F,pName,opt);
            end
        end
        
        function exportINP(F,filename)
            % Prepare voxel data for writing to INP file
            im = F.property.solid(1:(F.res(1)-1),1:(F.res(2)-1),1:(F.res(3)-1));
            [rows, cols, sli]  = size(im);
            ele_ind_vector = find(im == 1);
            num_ele = size(ele_ind_vector, 1);
            ele_temp = zeros(num_ele, 10);

            for j = 1: num_ele
                % subscript of element in im
                [ row, col, sli ] = ind2sub( [rows, cols, sli], ele_ind_vector(j) );

                % get linear index of eight corner of voxel(row,col,sli) in
                % nodecoor_list
                Lind_8corner = [
                    (col-1)*(rows+1) + row + (rows+1)*(cols+1)*(sli-1), ...
                    col*(rows+1) + row + (rows+1)*(cols+1)*(sli-1), ...
                    col*(rows+1) + row + 1 + (rows+1)*(cols+1)*(sli-1), ...
                    (col-1)*(rows+1) + row + 1 + (rows+1)*(cols+1)*(sli-1), ...
                    (col-1)*(rows+1) + row + (rows+1)*(cols+1)*sli, ...
                    col*(rows+1) + row + (rows+1)*(cols+1)*sli, ...
                    col*(rows+1) + row + 1 + (rows+1)*(cols+1)*sli, ...
                    (col-1)*(rows+1) + row + 1 + (rows+1)*(cols+1)*sli
                    ];
                % put Lind_8corner into
                ele_temp( j, : ) = [ ele_ind_vector(j), 1, Lind_8corner ];
            end
            ele_cell{1} = ele_temp;

            elements.E_ind = ele_temp(:,1);
            elements.E = ele_temp(:,3:10);

            % Generate unique nodes list
            node_ind_cell = cellfun( @(A) unique(A(:,3:10)), ele_cell, 'UniformOutput', 0 );
            unique_node_ind_v = unique( cell2mat( node_ind_cell ) );    % column vector

            % generate x y z coordinate of all nodes
            % can be accessed by X( row, col, sli ), Y( row, col, sli ),
            % Z( row, col, sli )
            [ X, Y, Z ] = ndgrid(F.xq, F.yq, F.zq);

            % reshape into vector
            % can be accessed by X(i), Y(i), Z(i)
            X = X(:);
            Y = Y(:);
            Z = Z(:);

            % extract certain nodes
            X = X( unique_node_ind_v );
            Y = Y( unique_node_ind_v );
            Z = Z( unique_node_ind_v );

            num_node = length( unique_node_ind_v );
            % temporary list for parfor
            temp_list = zeros( num_node, 3 );

            for i = 1: num_node
                temp_list( i, : ) = [ X(i), Y(i), Z(i) ];
            end
            nodeStruct.N_ind = unique_node_ind_v;
            nodeStruct.N = temp_list;
            elements.E_type = '*ELEMENT, TYPE=C3D8R, ELSET=TPMSDesigner-Elements';

            % Perform writing using template format
            writeINP(elements,nodeStruct,filename);
        end
    end
end

