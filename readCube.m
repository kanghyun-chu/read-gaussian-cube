function Cube = readCube( FileNames, flag )
% Kangh's readCube (2026.11.18) (Forked from Gary's repo.)
%   Reads Gaussian cube files into a struct array
%     Coordinates        Cube.axis_1 .axis_2 .axis_3
%     Spatial resolution Cube.resolution
%     Molecular config.  Cube.config
%     Volumetric data    Cube.cubedat
%     misc.              Cube.misc    (Reserved for extra data in Cube Header)
%   Convert the units from Bohr Radius to Angstroms upon demand
%
% In the case of non-orthogonal crystal axes, 
%  .resoultion is a 3 x 3 matrix composed of row-vectors: dAx1, dAx2, and dAx3.
%  .axis_123   is an array composed of [x;y;z]-coordinates in each column.


arguments
    FileNames string
    flag      double = 0
end

switch flag % Do unit conversion or not?
    case 0
        b0  = 1;
        b0cb= 1;
    case +1
        b0   = 0.52917721092;
        b0cb = b0^3;
    case -1
        b0   = 0.52917721092;
        b0cb = 1;
    otherwise
        b0 = 1; b0cb=1;
        warning("Invalid flag value. No spatial unit conversion is made.")
end

Cube = struct('filename',{}, ...
    'title'   ,{}, ...
    'subtitle',{}, ...
    'axis_1'  ,{}, ...
    'axis_2'  ,{}, ...
    'axis_3'  ,{}, ...
    'resolution',{}, ...
    'config'  ,{}, ...
    'cubedat' ,{}, ...
    'misc', {} );

for n=numel(FileNames):-1:1
    Cube(n).filename = FileNames(n); % Reverse order here for initializing struct object
end

% loop over FileNames
for n = 1:numel(FileNames)
    filename = FileNames(n);

    % Open file
    try
        if endsWith(filename, ".cube", "IgnoreCase",true)
            fid = fopen( filename );
        else
            fid = fopen( filename + ".cube" );
        end
    catch
        warning( "Failed to open file: %s", filename );
        continue
    end

    %% Reads headers lines 1-2
    Cube(n).title   = fgetl(fid);
    Cube(n).subtitle= fgetl(fid);

    %% read coordinate parameters Lines 3-6
    %   #_of_atoms , Origin_x , Origin_y, Origin_z
    %   #Grid_Axis1, dAxis1_x , dAxis1_y, dAxis1_z
    %   #Grid_Axis2, dAxis2_x , dAxis2_y, dAxis2_z
    %   #Grid_Axis3, dAxis3_x , dAxis3_y, dAxis3_z
    %
    % If #_of_atoms < 0, then extra line is collected and put to ".misc"
    %   For example, Orca generated .cube for Molecular orbital.

    param_coord = fscanf(fid, '%f %f %f %f \n', [4,4])'; % To have what-you-see-what-you-get, Do Transpose here.
    N_atoms = param_coord(1,1) ;
    Origin  = param_coord(1,2:4) * b0;
    NGr123  = param_coord(2:4,1)'; % Number of grid along 1,2,3 axis (in Row form)
    dAxis1  = param_coord(2,2:4) * b0; % Row vector
    dAxis2  = param_coord(3,2:4) * b0;
    dAxix3  = param_coord(4,2:4) * b0;

    dAx = [dAxis1;dAxis2;dAxix3]; % Each row is dAxis vector

    if isdiag(dAx) % Mutually orthogonal Axes
        %   Axis_{1,2,3} = {x,y,z}-coordinates
        dAx = reshape(diag(dAx), [1,3] );

        Cube(n).axis_1 = Origin(1) + (0:NGr123(1)-1)*dAx(1);
        Cube(n).axis_2 = Origin(2) + (0:NGr123(2)-1)*dAx(2);
        Cube(n).axis_3 = Origin(3) + (0:NGr123(3)-1)*dAx(3);
        Cube(n).resolution = dAx;

    else
        % .axis_# = [ Axis#_x(n=0)   ...   Axis#_x(n=end);
        %             Axis#_y(n=0)   ...   Axis#_y(n=end);
        %             Axis#_z(n=0)   ...   Axis#_z(n=end)]
        Cube(n).axis_1 = transpose(Origin + transpose((0:NGr123(1)-1)) * dAx(1,:));
        Cube(n).axis_2 = transpose(Origin + transpose((0:NGr123(2)-1)) * dAx(2,:));
        Cube(n).axis_3 = transpose(Origin + transpose((0:NGr123(3)-1)) * dAx(3,:));
        Cube(n).resolution = dAx;
    end

    % read molecule configuration lines 7-(7+number of molecules)

    if N_atoms >= 0
        config = fscanf(fid, '%f %f %f %f %f', [5,N_atoms])';
        config(:,3:5) = config(:,3:5) * b0; % Rescale atomic coordinates
        misc = [];
    else
        config = fscanf(fid, '%f %f %f %f %f', [5,abs(N_atoms)])';
        config(:,3:5) = config(:,3:5) * b0; % Rescale atomic coordinates
        misc = fscanf(fid, '%f %f', [2,1])'; % ORCA's Molecular orbital
    end
    Cube(n).config = array2table(config, ...
        'VariableNames',["Atomic #", "Charge","X","Y","Z"]);
    Cube(n).misc   = misc;

    %% Reads the volumetric data
    data = fscanf(fid, '%f' );

    if numel(data) ~= prod( NGr123 )
        warning(" Number of data point disagree in %s", filename)
        continue;
    end
    
    NGr_321 = NGr123(end:-1:1);
    data = reshape(data, NGr_321); % Array index ordering: axis_3,2,1
    
    Cube(n).cubedat = data / b0cb;

end
    Cube = reshape(Cube, size(FileNames) );
end
