function [matrix] = make_checks(SEGY,check)

% MAKE_CHECKERBOARD: Generates a checkerboard-pattern of perturbations 
%           applied to srModel fields u, anis_fractions, and aniso_phi.

% INPUT
%     srModel: Stingray structure
%       field: Model field in srModel to which the checkerboard will be
%                   applied.
%          xl: X-extent of checkerboard, [min-x, max-x]---km
%          yl: Y-extent of checkerboard, [min-y, max-y]---km
%          zl: Z-extent of checkerboard, [min-z, max-z]---km
%        cdim: Dimensions of checkers, [x-width, y-width, z-width]---km
%              Checker dimensions are rounded up if not an even multiple of
%              grid spacing.
%         mag: Magnitude of variation between adjacent checkers.
%              If 'field' is 'u' or 'anis_fraction', mag is the fractional
%              amount added or subtracted from blocks. Adjacent checkers will
%              have magnitudes of (1 + mag) and (1 - mag). The
%              specified srModel field is multiplied by the checker
%              magnitudes.
%
%              If 'field' is 'anis_phi', mag defines two orientations of
%              seismic anisotropy to prescribe to adjacent checkers.
%                   mag = [phi1, phi2]---radians
%
%
% OUTPUT
%     srModel: srModel with specified field modified by checkerboard
%      checks: checkerboard applied to srModel
%

% Unpack variables
xl = check.xl;
yl = check.yl;
zl = check.zl;
cdim = check.cdim;
mag = check.mag;

zl = fliplr(zl);

% Nodal dimensions of checks
cnx = ceil(cdim(1)/SEGY.gx);
cny = ceil(cdim(2)/SEGY.gy);
cnz = ceil(cdim(3)/SEGY.gz);

% Make ghead
SEGY.ghead = [SEGY.minx SEGY.miny SEGY.nx SEGY.ny SEGY.nz SEGY.gx SEGY.gy SEGY.gz];
% Indices of checkerboard extent
[io,jo,ko] = xyz2ijk(xl(1),yl(1),zl(1),SEGY.ghead);
[ie,je,ke] = xyz2ijk(xl(2),yl(2),zl(2),SEGY.ghead);

% Extend checkerboard limits such that an even number of blocks are used
% The checkerboard is extended in the direction of ii,ji,ki
ii = ie + cnx*ceil(length(io:ie)/cnx) - length(io:ie);
ji = je + cny*ceil(length(jo:je)/cny) - length(jo:je);
ki = ke + cnz*ceil(length(ko:ke)/cnz) - length(ko:ke);

% Create checkerboard foundation--4x4 alternating blocks in X,Y-plane
checks = cat(2,ones(cnx,cny),zeros(cnx,cny));
checks = cat(1,checks,~checks);

% Extend checkerboard in...
% X-direction,
checks = repmat(checks,ceil(0.5*length(io:ii)/cnx),1);
% Y-direction,
checks = repmat(checks,1,ceil(0.5*length(jo:ji)/cny));
% and Z-direction (need to invert pattern and extend two blocks deep)
checks = cat(3,repmat(checks,1,1,cnz),~repmat(checks,1,1,cnz));
checks = repmat(checks,1,1,ceil(0.5*(length(ko:ki))/cnz));
if mag < 1
    checks = (1+mag)*checks + (1-mag)*(~checks);
else
    checks = checks+mag + ~checks-mag;
end

xx = length(io:ie); dxx = xx*SEGY.gx;
yy = length(jo:je); dyy = yy*SEGY.gy;
zz = length(ko:ke)-1; dzz = zz*SEGY.gz;
if round(dxx)/SEGY.gx ~= xx
    while round(dxx) > xx*SEGY.gx
        ie = ie+1;
        xx = length(io:ie);
    end
    while round(dxx) < xx*SEGY.gx
        ie = ie-1;
        xx = length(io:ie);
    end
    
end

if round(dyy)/SEGY.gy ~= yy 
    while round(dyy) > yy*SEGY.gy
        je = je+1;
        yy = length(jo:je);
    end
    while round(dyy) < yy*SEGY.gy
        je = je-1;
        yy = length(jo:je);
    end
end

if (floor(dzz*10)/10)/SEGY.gz ~= zz
    while roundn(dzz,-1) > zz*SEGY.gz
        ke = ke+1;
        zz = length(ko:ke);
    end
    while round(dzz) < zz*SEGY.gz
        ke = ke-1;
        zz = length(ko:ke);
    end
end


checks = checks(1:length(io:ie),(1:length(jo:je)-1),1:length(ko:ke)-1);

matrix = ones(SEGY.nx, SEGY.ny, SEGY.nz);
matrix(io:ie,jo:je-1,ko:ke-1) = checks;


    
end
