from __future__ import print_function, absolute_import
from scipy.interpolate import interp1d
import underworld as uw
#from UWGeodynamics import nd
from underworld import UWGeodynamics as GEO
from mpi4py import MPI as _MPI

comm = _MPI.COMM_WORLD
rank = comm.rank
size = comm.size

class FreeSurfaceProcessor(object):
    """FreeSurfaceProcessor"""

    def __init__(self, model):
        """Create a Freesurface processor

        Parameters
        ----------

        model : UWGeodynamics Model

        """
        self.model = model

        minCoord = self.model.fsminCoord
        maxCoord = self.model.fsmaxCoord

        # Initialize model mesh
        self._init_mesh = uw.mesh.FeMesh_Cartesian(elementType=self.model.elementType,
                                                   elementRes=(self.model.fsmeshnx,self.model.fsmeshny),
                                                   minCoord=minCoord,
                                                   maxCoord=maxCoord,
                                                   periodic=self.model.periodic)
                                     
        # Create the tools
        self.TField = self._init_mesh.add_variable(nodeDofCount=1)
        self.TField.data[:, 0] = self._init_mesh.data[:, 1].copy()

        self.top = self._init_mesh.specialSets['Top_VertexSet'] 
        self.bottom = self._init_mesh.specialSets['Bottom_VertexSet'] 

        # Create boundary condition
        self._conditions = uw.conditions.DirichletCondition(
            variable=self.TField,
            indexSetsPerDof=(self.top + self.bottom,))

        # Create Eq System
        self._system = uw.systems.SteadyStateHeat(
            temperatureField=self.TField,
            fn_diffusivity=1.0,
            conditions=self._conditions)

        self._solver = uw.systems.Solver(self._system)

    def _solve_sle(self):
        self._solver.solve()

    def _advect_surface(self, dt):

        if self.top:
            # Extract top surface
            x = self.model.mesh.data[self.model.top_wall.data][:, 0]
            y = self.model.mesh.data[self.model.top_wall.data][:, 1]

            # Extract velocities from top
            vx = self.model.velocityField.data[self.model.top_wall.data][:, 0]
            vy = self.model.velocityField.data[self.model.top_wall.data][:, 1]

            # Advect top surface
            x2 = x + vx * GEO.nd(dt)
            y2 = y + vy * GEO.nd(dt)

            # Spline top surface
            f = interp1d(x2, y2, kind='cubic', fill_value='extrapolate')

            self.TField.data[self.top.data, 0] = f(x)
            
        comm.Barrier()
        self.TField.syncronise()
        
    def _update_mesh(self):
        with self.model.mesh.deform_mesh():
        	self.model.mesh.data[(self.model.fsmeshnx+1)*(self.model.fsmeshny1):, -1] = self.TField.data[:, 0].copy()

    def solve(self, dtime):
        """ Advect free surface through dt and update the mesh """

        # First we advect the surface
        self._advect_surface(dtime)
        # Then we solve the system of linear equation
        self._solve_sle()
        # Finally we update the mesh
        self._update_mesh()
