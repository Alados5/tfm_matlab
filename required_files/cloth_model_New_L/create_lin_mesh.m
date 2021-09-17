function MeshPos = create_lin_mesh(lSide, nSide, cpt, angle)
% Creates a cloth mesh in space. Each row of the output is a node,
% with columns corresponding to coordinates [x,y,z]
% The cloth center is cpt, and its side length is lSide, which can be
% a number (square cloth) or a vector of two (rectangular)
% The mesh size is set by nSide, which can also be a number or two
% By default, the mesh is created on the XZ plane, but can be rotated
% along the Z axis with angle
%
% - Original Author: Adrià Luque, adria.alados@gmail.com
% - Last review: September 2021


lx = lSide(1);
if length(lSide) == 1
    lz = lx;
else
    lz = lSide(2);
end

nx = nSide(1);
if length(nSide) == 1
    nz = nx;
else
    nz = nSide(2);
end

[px,py,pz] = meshgrid(linspace(-lx/2,lx/2,nx), ...
                      0, ...
                      linspace(-lz/2,lz/2,nz));

px=permute(px,[2,3,1]);
py=permute(py,[2,3,1]);
pz=permute(pz,[2,3,1]);
px=reshape(px,[nx*nz,1]);
py=reshape(py,[nx*nz,1]);
pz=reshape(pz,[nx*nz,1]);

MeshPosXZ = [px py pz];

RotM = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];
MeshPos = (RotM * MeshPosXZ')' + cpt;

end