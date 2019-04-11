''' Python script to plot the hydro variables '''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

outdir = "cyl_Re80.out/"

# Output all frames
frames = range(10)
# Output selected frames
# frames = [0,10,20,30]

# Load one .txt file and return the absolute value of velocity and density
def load_file(fr):
	filename = "fr_%04d"%fr
	contents = np.loadtxt(outdir+filename+".txt")
	header = contents[0,:]
	nxb = int(header[1])
	nyb = int(header[2])
	rho = contents[1:nxb*nyb+1,0].reshape(nyb,nxb)
	ux = contents[1:nxb*nyb+1,1].reshape(nyb,nxb)
	uy = contents[1:nxb*nyb+1,2].reshape(nyb,nxb)
	usq = np.sqrt(ux*ux+uy*uy)
	X = np.arange(0,nyb-2,1)
	Y = np.arange(0,nxb-2,1)
	X, Y = np.meshgrid(X, Y)
	return rho[1:nyb-1,1:nxb-1],usq[1:nyb-1,1:nxb-1],X,Y

# Get limits for colorbar range
def get_lims(rho,usq):
	bar1 = np.zeros_like(rho)
	bar1[rho==0] = 1
	bar1[bar1==0] = np.nan
	bar2 = np.zeros_like(usq)
	bar2[usq==0] = 1
	bar2[bar2==0] = np.nan
	usqmin = np.min(usq)
	usqmax = np.max(usq)
	rhomin = np.min(rho[rho>0])
	rhomax = np.max(rho)
	return bar1,bar2,usqmin,usqmax,rhomin,rhomax

# Open the last frame to get colorbar limits
rho,usq = load_file(frames[-1]-1)[0:2]
bar1,bar2,usqmin,usqmax,rhomin,rhomax = get_lims(rho,usq)

for fr in frames:
	filename = "fr_%04d"%fr
	rho,usq,X,Y = load_file(fr)
	fig, ax = plt.subplots(figsize=(15,6))
	# print("Frame: ",fr)
	# print(usq[1:ny-1,1:nx-1])
	cs = ax.imshow(usq,cmap=cm.magma)
	ax.imshow(bar2,cmap=cm.nipy_spectral_r)
	cs.set_clim(0,np.round(usqmax,3))
	cbar = fig.colorbar(cs, ticks=[0,np.round(usqmax/2,3),np.round(usqmax,3)], orientation='horizontal')
	cbar.ax.tick_params(labelsize=28)
	ax.set_xlabel(r'$x$', fontsize=16)
	ax.set_ylabel(r'$y$', fontsize=16)
	ax.set_title('Velocity at Frame %04d' % fr, fontsize=20)
	plt.savefig(outdir+'%04d.png' % fr)
	plt.close('all')