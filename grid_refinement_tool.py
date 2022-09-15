import netCDF4

import numpy

from meshkernel import MeshKernel, GeometryList, Mesh2dLocation,AveragingMethod

from meshkernel import RefinementType, MeshRefinementParameters

from matplotlib import pyplot, collections

def crop_mesh( m, polygon):
        
    outside=m.mesh2d_get_nodes_in_polygons(polygon, False)
  
    for n in outside:
      m.mesh2d_delete_node(n)
  
    m.mesh2d_get_hanging_edges() # needed otherwise deletes all edges
    m.mesh2d_delete_hanging_edges()

def refine_mesh( m, polygon, max_refinement_level=1):
         
    parameters=MeshRefinementParameters(refine_intersected=False, 
                                        use_mass_center_when_refining=False,
                                        min_face_size=1.,
                                        refinement_type=1,#RefinementType.REFINEMENT_LEVELS,
                                        connect_hanging_nodes=False,
                                        account_for_samples_outside_face=False, max_refinement_iterations=max_refinement_level)
             
    m.mesh2d_refine_based_on_polygon(polygon, parameters)

def clone_export_dataset(filename, d, NetNode_x, NetNode_y, NetNode_z, NetLink, NetLinkType):
    """ export mesh to netcdf dataset using d as template """
    
    ds=netCDF4.Dataset(filename, "w", format="NETCDF3_CLASSIC")
    
    # copy attributes
    for name in d.ncattrs():
      ds.setncattr(name, d.getncattr(name))
    
    # maybe add better description to history
    ds.setncattr("history", "cutout;"+d.getncattr("history"))
    
    # create dimensions
    ds.createDimension("nNetNode", size= len(NetNode_x))
    ds.createDimension("nNetLink", size= len(NetLink))
    ds.createDimension("nNetLinkPts", size=2)
    #~ ds.createDimension("nNetElem", size=?)
    
    for name, variable in d.variables.items():
      print(name, variable.dimensions)
      var=ds.createVariable(name, variable.datatype, variable.dimensions)
      for aname in variable.ncattrs():
        var.setncattr(aname, variable.getncattr(aname))
    
    ds["NetNode_x"][:]=NetNode_x
    ds["NetNode_y"][:]=NetNode_y
    ds["NetNode_z"][:]=NetNode_z
    ds["NetLink"][:,:]=NetLink+1 # back to 1 based indexing 
    ds["NetLinkType"][:]=NetLinkType
        
    ds.close()

# data structure to contain mesh information of netfile
# including a reference to the original dataset
class Grid:
    def __init__(self, NetNode_x, NetNode_y, NetNode_z, NetLink, NetLinkType, parentds):
        self.NetNode_x=NetNode_x
        self.NetNode_y=NetNode_y
        self.NetNode_z=NetNode_z
        self.NetLink=NetLink
        self.NetLinkType=NetLinkType
        self.ds=parentds
        
    def save_to_file(self, filename):
        clone_export_dataset(filename, self.ds, self.NetNode_x, self.NetNode_y, self.NetNode_z,
           self.NetLink, self.NetLinkType)

    def plot(self,filename="test.png"):
        fig, ax = pyplot.subplots()
      
        x=self.NetNode_x
        y=self.NetNode_y
        z=self.NetNode_z
      
        pyplot.scatter(x,y,c=z, cmap="jet", s=1.)  
      
        pyplot.savefig(filename)
      
# so what happens around the wrappingline?
class GridRefinementTool:

    def __init__(self, filename,spherical_coordinates=True, interpolator=None):
        self.d=netCDF4.Dataset(filename)
        # probably can be derived from d
        self.spherical_coordinates=spherical_coordinates
  
        self.read_data()

        if interpolator is None:
            interpolator=self.triangulation_interpolator
        self.interpolator=interpolator

    def read_data(self):
        self.NetNode_x=self.d["NetNode_x"][:]
        self.NetNode_y=self.d["NetNode_y"][:]
        self.NetNode_z=self.d["NetNode_z"][:]
        self.NetLink=self.d["NetLink"][:]-1 # zero based indexing
        self.NetLinkType=self.d["NetLinkType"][:]

        assert(numpy.all(self.NetLinkType==2))

    def import_net_dataset_to_meshkernel(self):
                
        m=MeshKernel(is_geographic=self.spherical_coordinates)
        
        _nodes=zip(self.NetNode_x,self.NetNode_y)
        _edges=self.NetLink
        
        nodes=[]
        for x in _nodes:
          n=m.mesh2d_insert_node(*x)
          nodes.append(n)
          
        edges=[]
        for x in _edges:
          e=m.mesh2d_insert_edge(*x)
          edges.append(e)
                
        return m 
        
    @property
    def grid(self):
        return Grid(self.NetNode_x, self.NetNode_y, self.NetNode_z, self.links, self.ltype, self.d)

    def refine_around(self, location, cutout_size, refinement_size, refinement_level=3):
        
        x1=[location[0]-cutout_size/2, location[0]+cutout_size/2]
        y1=[location[1]-cutout_size/2, location[1]+cutout_size/2]

        x2=[location[0]-refinement_size/2, location[0]+refinement_size/2]
        y2=[location[1]-refinement_size/2, location[1]+refinement_size/2]
        
        
        polygon1=GeometryList(numpy.array([x1[0],x1[1],x1[1],x1[0],x1[0]]),
                       numpy.array([y1[0],y1[0],y1[1],y1[1],y1[0]]))
        
        polygon2=GeometryList(numpy.array([x2[0],x2[1],x2[1],x2[0],x2[0]]),
                       numpy.array([y2[0],y2[0],y2[1],y2[1],y2[0]]))
                
        m=self.import_net_dataset_to_meshkernel()
        
        crop_mesh(m, polygon1)
        
        if len(m.mesh2d_get().node_x)==0:
            return None
        
        if refinement_level>0:
            refine_mesh(m, polygon2, refinement_level)


        result=self.interpolator(m)
  
        NetNode_x=result.x_coordinates
        NetNode_y=result.y_coordinates
        NetNode_z=result.values
      
        mesh2d=m.mesh2d_get()
      
        links=mesh2d.edge_nodes.reshape((-1,2))
      
        ltype=numpy.zeros(len(links))+2 # hardcodes

        return Grid(NetNode_x, NetNode_y,NetNode_z, links, ltype, self.d)

    def triangulation_interpolator(self, m):
        
        geometrylist=GeometryList(self.NetNode_x, self.NetNode_y, self.NetNode_z)

        # returns Geometry list with interpolated values on mesh locations
        return m.mesh2d_triangulation_interpolation(geometrylist, Mesh2dLocation(1))

        

if __name__=="__main__":
  
    #~ location=[(145.224609375+152.783203125)/2, (-22.8515625+ -15.0)/2 ]
    #~ cutout_size=(152.783203125-145.224609375)

    location=[4.,52.]
    cutout_size=14

    refinement_size=cutout_size-3.
    refinement_level=0

    #~ g=GridRefinementTool("test.nc")
    g=GridRefinementTool("step11_global_net.nc")

    refinement=g.refine_around( location, cutout_size, refinement_size, refinement_level)

    if refinement is None:
        print("empty region selected")
    else:
        refinement.save_to_file("test_refinement.nc")
        refinement.plot("refinement.png")
