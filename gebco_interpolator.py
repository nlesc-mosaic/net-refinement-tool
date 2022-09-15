import netCDF4
import numpy

from meshkernel import GeometryList

GEBCO_DATAFILE="/home/inti/code/mosaic/gebco_data/GEBCO_2020.nc"

data=netCDF4.Dataset(GEBCO_DATAFILE)

Nlon=data["lon"].shape[0]
Nlat=data["lat"].shape[0]

dlon=360./Nlon
dlat=180./Nlat

lon0=data["lon"][0]
lat0=data["lat"][0]

elev=data["elevation"]

assert dlon==dlat

class gebco_interpolator:
    def __init__(self, resolution=None):
        # eventually resolution indicates resolution to which is averaged
        resolution_factor=1
        if resolution_factor==1:
            self.elev=data["elevation"]
            self.Nlon=Nlon
            self.Nlat=Nlat
            self.lon0=lon0
            self.lat0=lat0
            self.dlon=dlon
            self.dlat=dlat

    def interpolate(self, lon,lat):
        """
          interpolate GEBCO data linearly. returns elev data by bilinear interpolation
          of the input lon, lat points
        """
        
        lon=numpy.mod(lon+180, 360.)-180
        if numpy.any(lat<-90) or numpy.any(lat>90):
            raise Exception("invalid lattitude in input")
    
        ix=(lon-self.lon0)/self.dlon
        wx=1-(ix-numpy.floor(ix))
        ix=numpy.floor(ix).astype(int)
        ix1=ix+1
    
        ix[ix==-1]=self.Nlon-1
        ix1[ix==self.Nlon]=0
    
        iy=(lat-self.lat0)/self.dlat
        wy=1-(iy-numpy.floor(iy))
        iy=numpy.floor(iy).astype(int)
        iy1=iy+1
        
        iy=numpy.clip(iy,0, self.Nlat)
        iy1=numpy.clip(iy1,0, self.Nlat)
                
        result=         wx*wy*self.elev[iy ,ix ] + \
                   (1.-wx)*wy*self.elev[iy ,ix1] + \
                   wx*(1.-wy)*self.elev[iy1,ix ] + \
              (1.-wx)*(1.-wy)*self.elev[iy1,ix1]
    
        return result
    
    def interpolator(self, m):
        mm=m.mesh2d_get()
        x=mm.node_x
        y=mm.node_y
        values=self.interpolate(x,y)
        return GeometryList(x,y,values)
  
if __name__=="__main__":
    lon=numpy.array([-180+dlon/2+dlon/10])
    lat=numpy.array([-90+dlat/2+dlon/5])
    interpolate(lon,lat)
