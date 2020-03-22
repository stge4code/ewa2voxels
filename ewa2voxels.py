import struct
import sys,os
import numpy
import argparse
from pyevtk.hl import pointsToVTK 


class ReadEOFException(Exception):
    def __init__(self):
        pass

class Ewald3D():
    def __init__(self,ewafilename):
        self.ewafile = open(ewafilename, 'rb')  
        self.vtufilename = os.path.splitext(ewafilename)[0]
        self.datfilename = os.path.splitext(ewafilename)[0] + ".dat"
        self.FileVersion = None
        self.NumberOfRuns = None
        self.OriginPosition = None
        self.Scale = None
        self.VoxelAdditionalDataSize = None
        self.Ewald3DEmergencyData1 = None
        self.Ewald3DEmergencyData2 =  None
        self.IsMinMaxDataPresent = None
        self.Ewald3DEmergencyData3 = None
        self.Ewald3DEmergencyData4 = None
        self.MinIntensity = None
        self.MaxIntensity = None
        self.AverageIntensity = None
        self.SigmaOfIntensity = None 
        self.MaskSize = None
        self.VoxelsData = None
        self.VoxelsDataCounter = None

    def read(self,iformat,isize):
        raw_value = self.ewafile.read(isize)        
        if isize != len(raw_value):
            raise ReadEOFException
        value = struct.unpack(iformat, raw_value)
        return value
  
    def read_header(self):
        EEDSU = 50
        self.FileVersion = self.read('1L',4)[0]
        self.NumberOfRuns =self.read('1L',4)[0]
        self.OriginPosition = self.read('3L',3*4)
        self.Scale = self.read('<1f',4)[0]
        self.VoxelAdditionalDataSize = self.read('<1L',4)[0]
        self.Ewald3DEmergencyData1 = self.read('<{0:.0f}B'.format(EEDSU),EEDSU)
        self.Ewald3DEmergencyData2 =  self.read('<{0:.0f}h'.format(EEDSU-1),(EEDSU-1)*2)
        self.IsMinMaxDataPresent = self.read('<2h',4)
        self.Ewald3DEmergencyData3 = self.read('<{0:.0f}L'.format(EEDSU),(EEDSU)*4)
        self.Ewald3DEmergencyData4 = self.read('<{0:.0f}d'.format(EEDSU-2),(EEDSU-2)*8)
        self.MinIntensity = self.read('1f',4)[0]
        self.MaxIntensity = self.read('1f',4)[0]
        self.AverageIntensity = self.read('1f',4)[0]
        self.SigmaOfIntensity = self.read('1f',4)[0]
        self.MaskSize = int(numpy.ceil(self.NumberOfRuns/8))

    def read_voxels(self):        
        VoxelsData = []
        VoxelsDataCounter = 0 
        while True:
            try:
                vc = self.read('1L',4)[0]
                bvc = "{0:b}".format(vc)
                VoxelCoordinates = [int(bvc[-10:],2), int(bvc[-20:-10],2),int(bvc[-30:-20],2)]          
                IntensityValue = self.read('1f',4)[0]
                if self.VoxelAdditionalDataSize:
                    self.read('<{0:.0f}B'.format(self.VoxelAdditionalDataSize),self.VoxelAdditionalDataSize)[0]
                if self.MaskSize: self.read('<{0:.0f}B'.format(self.MaskSize),self.MaskSize)[0]
                VoxelsData.append([*VoxelCoordinates,IntensityValue])            
                VoxelsDataCounter += 1
            except ReadEOFException:
                break  
        self.VoxelsData = numpy.array(VoxelsData)
        self.VoxelsDataCounter = VoxelsDataCounter
        return  self.VoxelsData

    def print_header(self):
        print("Header Data:") 
        print(" File vesrsion: {0:.0f}".format(self.FileVersion)) 
        print(" Origin Position: {0:.0f} {1:.0f} {2:.0f}".format(*self.OriginPosition)) 
        print(" Scale: {0:.3f}".format(self.Scale)) 
        print(" Number Of Runs: {0:.0f}".format(self.NumberOfRuns)) 
        print(" Min Intensity: {0:.3f}".format(self.MinIntensity)) 
        print(" Max Intensity: {0:.3f}".format(self.MaxIntensity)) 
        print(" Average Intensity: {0:.3f}".format(self.AverageIntensity)) 
        print(" Sigma Of Intensity: {0:.3f}".format(self.SigmaOfIntensity)) 
 
    def print_voxels(self):
        print("Voxels Data:")         
        print(" {0:.0f} <= X <= {1:.0f}".format(numpy.min(self.VoxelsData[:,0]),numpy.max(self.VoxelsData[:,0]))) 
        print(" {0:.0f} <= Y <= {1:.0f}".format(numpy.min(self.VoxelsData[:,1]),numpy.max(self.VoxelsData[:,1])))
        print(" {0:.0f} <= Z <= {1:.0f}".format(numpy.min(self.VoxelsData[:,2]),numpy.max(self.VoxelsData[:,2])))
        print(" Voxels Count: {0:.0f}".format(self.VoxelsDataCounter)) 

 
    def export_vtu(self):     
        print("Creating VTU data file...",end="\r")   
        pointsToVTK(self.vtufilename, 
            numpy.ascontiguousarray(self.VoxelsData[:,0]),
            numpy.ascontiguousarray(self.VoxelsData[:,1]), 
            numpy.ascontiguousarray(self.VoxelsData[:,2]), 
            data = {"Intensity" : numpy.ascontiguousarray(self.VoxelsData[:,3])})
        print("VTU data file: " + self.vtufilename + ".vtu") 
    
   
    def export_xyzI(self):      
        print("Creating DAT file...",end="\r")   
        with open(self.datfilename,"w") as datfile:
            for item in self.VoxelsData:
                datfile.write("{0:4.0f} {1:4.0f} {2:4.0f} {3:.5f}\n".format(*item)) 
        print("(X,Y,Z,I) data file: " + self.datfilename) 



    def __del__(self):
        self.ewafile.close()  

            
if __name__ == "__main__": 
    parser = argparse.ArgumentParser(description='Convert .ewa file to *.vtu and (X,Y,Z,I) array.')
    parser.add_argument('inputewafile', type=str, help='*.ewa file name')
    args = parser.parse_args()
    ewa = Ewald3D(args.inputewafile)
    ewa.read_header()
    ewa.read_voxels()
    ewa.print_header()
    ewa.print_voxels()
    ewa.export_vtu()
    ewa.export_xyzI()

    
