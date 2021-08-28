function MeshPos = create_lin_mesh(lSide, nSide, cpt, angle)

nx = nSide(1);
if length(nSide) == 1
    nz = nx;
else
    nz = nSide(2);
end

[px,py,pz] = meshgrid(linspace(-lSide/2,lSide/2,nx), ...
                      0, ...
                      linspace(-lSide/2,lSide/2,nz));

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