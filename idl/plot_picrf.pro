pro plot_picrf

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

!p.multi = [0,2,2]
!p.charSize = 2.0

plot, x[0,*], $
		yRange = [0,1], $
		psym = 2

oPlot, x[1,*], $
		psym = 4

plot, pe
oplot, ke

plot, ke + pe 

surface, rho, zRange = [-1, 1], zStyle = 0, charSize = 2.0, font=0

stop
for i=0,n_elements(x[0,*])-1 do begin

!p.multi = [0,1,2]
plot, x[*,i], vx[*,i], $
		psym = 1, $
		yRange = [-1e4,1e4], $
		yStyle = 1

!p.multi = [2,2,2]
vx_hist = histogram ( vx[*,i], nBins = 40, locations = locs )
plot, locs, vx_hist, $
		yRange = [0,250], $
		xRange = [-1e4,1e4], $
		xStyle = 1, $
		yStyle = 1

plot, pe[0:i], yRange = [0,1e6]
oPlot, ke[0:i]
oPlot, pe[0:i] + ke[0:i]

wait, 0.005

endfor

stop
end
