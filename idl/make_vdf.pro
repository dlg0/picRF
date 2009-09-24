pro make_vdf

fileName	= 'output/output.nc'
cdfId	= ncdf_open ( fileName, /noWrite )
ncdf_varGet, cdfId, 'x', x 
ncdf_varGet, cdfId, 'vx', vx 
ncdf_varGet, cdfId, 'rho', rho
ncdf_varGet, cdfId, 'phi', phi
ncdf_varGet, cdfId, 'Ex', Ex 
ncdf_varGet, cdfId, 'dx', dx
ncdf_varGet, cdfId, 'weight', weight
ncdf_close, cdfId

e   = 1.60217646e-19
me  = 9.10938188e-31
mi  = 1.67262158e-27
e0  = 8.85418782e-12

ke	= 0.5 * me * total ( vx^2 * weight, 1 ) / e * 1e-3
pe	= dx / 2.0 * total ( rho * phi, 1 ) / e * 1e-3
energy = ke + pe



stop
end
