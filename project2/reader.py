import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

reader = vtk.vtkXMLUnstructuredGridReader()
# reader.SetFileName("BUFFY/RBC_0_Re10/DIRDATA_BloodFlow/VTK/T0000232000.pvtu")
reader.SetFileName("data/RBC_0_Re10/DIRDATA_BloodFlow/VTK/T0000100000.pvtu")
reader.Update()
data = reader.GetOutput()
points = data.GetPoints()
x = vtk_to_numpy(points.GetData())
print(x)
