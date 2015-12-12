function [meshNodes,meshElems] = GenerateDogBoneMesh(w,h,L,r,meshSize)

addpath(genpath('distmesh'));
addpath(genpath('postprocessing'));

fd=@(p) ddiff(ddiff(ddiff(ddiff(dunion(dunion(dunion(dunion(drectangle(p,-w,w,0,h),drectangle(p,-w+r,w-r,h,h+r)),drectangle(p,-w+2*r,w-2*r,h,h+L)),drectangle(p,-w+r,w-r,h+L-r,h+L)),drectangle(p,-w,w,h+L,2*h+L)),dcircle(p,-w+r,h+r,r)),dcircle(p,w-r,h+r,r)),dcircle(p,w-r,h+L-r,r)),dcircle(p,-w+r,h+L-r,r));

[meshNodes,meshElems]=distmesh2d(fd,@huniform,meshSize,[-2*w,-w;2*w,2*L],[-w,0;-w,h;-w,2*h+L;-w,h+L;w,0;w,h;w,2*h+L;w,h+L]);

end

