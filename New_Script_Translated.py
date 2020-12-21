@ File    (label = "Input directory", style = "directory") srcFile
@ File    (label = "Output directory", style = "directory") dstFile
@ String  (label = "File extension", value=".oir") ext
@ String  (label = "File name contains", value = "") containString
@ boolean (label = "Keep directory structure when saving", value = true) keepDirectories
@ Integer (label = "Maximum size", value=100) cellradius #maximum expected cell radius in µm
@ Integer (label = "nucleus channel", value=2) x #nucleus channel
@ Integer (label = "brightfield channel", value=3) y #brightfield channel

import os
import ij
from ij import IJ, ImagePlus
import os.path as path
from ij import WindowManager as WM  
from ij.measure import ResultsTable
from ij.text import TextWindow
import math
from ij.plugin import ImageCalculator as ic
import ij.process.ImageProcessor as ip
from ij.gui import HistogramWindow
from ij.plugin import RGBStackMerge
from ij.process import ImageProcessor
from ij.plugin.filter import MaximumFinder
from jarray import array
from ij import Prefs
from ij.plugin.filter import Binary
from ij.plugin.frame import RoiManager
from ij.gui import Roi
from ij.plugin.filter import ParticleAnalyzer
#from ij.gui import PolygonRoi

f = "[Cell_summaries]" #title of the produced 
seriesnum = 1 #Number of series in the leica image file
counter = 1


IJ.run("Close All");
#setBatchMode(true); #prevents image windows from opening while the script is running
print("\\Clear");
print("\\Update:<Processing started>")
IJ.run("Clear Results") 

def processFolder():
  srcDir = srcFile.getAbsolutePath()
  dstDir = dstFile.getAbsolutePath()
  for root, directories, filenames in os.walk(srcDir):
    filenames.sort();
    for filename in filenames:
      # Check for file extension
      if not filename.endswith(ext):
        continue
      # Check for file name pattern
      if containString not in filename:
        continue
      cellSegmentation(srcDir, dstDir, root, filename, keepDirectories)

  
def save_all(srcDir, dstDir, filename, localfile, keepDirectories, imageID): 
  saveDir = currentDir.replace(srcDir, dstDir) if keepDirectories else dstDir
  if not os.path.exists(saveDir):
    os.makedirs(saveDir)
  print "Saving to", saveDir
  IJ.selectWindow("Results for PA ")
  IJ.saveAs("Text", ""+saveDir+"\\"+localfile+".xls")
  IJ.selectWindow(imageID);
  imp=IJ.getImage()
  IJ.run("Duplicate...", "duplicate channels="+str(y)+"") #brightfield channel
  IJ.run("Z Project...", "projection=[Sum Slices]")
  rm= RoiManager(True)
  rm.runCommand("Show All without labels")
  rm.runCommand("Set Color", "red")
  rm.runCommand("Set Line Width", str(4))
  imp=IJ.getImage()
  IJ.saveAs(imp,"Tiff", ""+saveDir+"\\"+localfile+"_OV.tif")
  rm.runCommand("Show None")
  IJ.saveAs(imp,"Tiff", ""+saveDir+"\\"+localfile+".tif")
  IJ.selectWindow("IMAGE")
  imp=IJ.getImage()
  IJ.saveAs(imp, "Tiff", os.path.join(saveDir, filename))
  ROInumber = rm.getCount()
  if(ROInumber>0):
	RoiManager("Delete")

def cellSegmentation(srcDir, dstDir, currentDir, filename, keepDirectories):
  print "Processing:"  
  # Opening the image
  print "Open image file", filename
  imp = IJ.openImage(os.path.join(currentDir, dstDir))
  # Put your processing commands here!
  localinput=srcDir.replace("/", "\\")
  saveDir = localinput.replace(srcDir, dstDir)
  string="."
  dotIndex=filename.find(string)
  localfile= filename[0:dotIndex]
  print(localfile)
  IJ.run("New... ", "name="+f+" type=Table")
  print(f,"\\Headings:Cell\tarea\tCirc\tAR\tRoundness\tMaximum")
  IJ.run("Bio-Formats", "open=[" + localinput + os.path.sep + filename +"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT")
  IJ.open()
  idd= WM.getIDList();
  imageID= idd[0];
  IJ.run("Clear Results")
  WM.getImage(imageID)
  IJ.run("Duplicate...", "duplicate channels="+str(x)+"") #Nucleus channel #took away x
  IJ.run("Z Project...", "projection=[Standard Deviation]");#picture for frame detection
  IJ.run("8-bit");
  IJ.run("Duplicate...", "title=IMAGE");#frame
  IJ.run("Duplicate...", "title=SUBTRACT");#Background subtraction mask (for frame and watershed)
  imp=IJ.getImage()
  pixelWidth=imp.getWidth()
  pixelWidth=pixelWidth/1647.89
  pixelHeight= imp.getHeight()
#create subtraction mask, applying constraining maximum (step I)
  IJ.selectWindow("SUBTRACT")
  nResults=imp.getStatistics()
  row = nResults
  rt_exist = WM.getWindow("Results")
  if rt_exist==None:
    rt= ResultsTable()
  else:
    rt = rt_exist.getTextPanel().getOrCreateResultsTable()
  rt.setValue("Max ", 0, row.max) #text file
  rt.show("Results")
  u=math.floor(row.mean*3)
  IJ.run("Max...","value="+str(u)) #constraining maximum of 3-fold mean to reduce effect of extreme values during subtraction
		#gaussian blurring (step II)
  IJ.run("Gaussian Blur...", "sigma=100 scaled") #blurring for subtraction mask

  IJ.selectWindow("IMAGE")
  pxrollrad = cellradius/pixelWidth; #rolling ball radius in pixels needed (= predefined cell radius[µm]/pixelsize[µm/px])
  IJ.run("Subtract Background...", "rolling="+str(pxrollrad)+"")
  IJ.run("Gaussian Blur...", "sigma=2 scaled") #reduces punctate character of grayscale image '
  IM=IJ.selectWindow("IMAGE")
  SUB=IJ.selectWindow("SUBTRACT")
  ic().run("SUBTRACT", IM, SUB) #just subtracts two images
  IJ.selectWindow("IMAGE") #see how to call
  IJ.run("Duplicate...", "title=AND")#watershed
  IJ.run("Duplicate...", "title=CHECK")#for checking if maxima exist within selection later
  
#Apply threshold to get binary image of cell borders (step IV)
  IJ.selectWindow("IMAGE")
  imp = IJ.getImage()  # the current image
  imp.getProcessor().setThreshold(1, 255, ImageProcessor.NO_LUT_UPDATE)
  IJ.run("Subtract Background...","...")
  IJ.run("Convert to Mask", "method=Default background=Dark only black")
  IJ.run("Fill Holes")

#Create watershed line image (step V)
  IJ.selectWindow("AND")
  IJ.run("Gaussian Blur...", "sigma=2 scaled")
  imp=IJ.getImage()
  pixelWidth=imp.getWidth()
  pixelWidth=pixelWidth/1647.89
  pixelHeight= imp.getHeight()
  # Saving the image
  nResults=imp.getStatistics()
  row = nResults
  rt.setValue("Max ", 1, row.max) #text file
  nBins = 256
  Hist = HistogramWindow("Histogram",imp,nBins)
  Table = Hist.getResultsTable()
  Counts = Table.getColumn(1)
  #mean gray value of pixels belonging to cells needed (i.e. mean of ONLY non-zero pixel)
  Sum = 0 #all counts
  CV = 0 #weighed counts (= counts * intensity)
  for i in range(0, len(Counts)): #starting with 1 instead of 0. -> 0 intensity values are not considered.
    Sum += Counts[i]
    CV += Counts[i]*i
  m = (CV/Sum)
  m=math.floor(m)
  l = math.floor(2*m) #Maxima need to be at least twice the intensity of cellular mean intensity
  IJ.run("Find Maxima...", "noise="+str(l)+" output=[Segmented Particles] exclude") #watershedding

#Combine watershed lines and cell frame (step VI) 
  IJ.selectWindow("IMAGE")
  imp=IJ.getImage()
  imp.getProcessor().setThreshold(1, 255, ImageProcessor.NO_LUT_UPDATE)
  IJ.run(imp, "Watershed", "") #useful
  imp = IJ.getImage()
  ip = imp.getProcessor()
  segip = MaximumFinder().findMaxima( ip, 1, ImageProcessor.NO_THRESHOLD, MaximumFinder.SEGMENTED , False, False)
  segip.invert()
  segimp = ImagePlus("seg", segip)
  segimp.show()
  mergeimp = RGBStackMerge.mergeChannels(array([segimp, None, None, imp, None, None, None], ImagePlus), True)
  mergeimp.show()
  pa_exist = WM.getWindow("Results for PA")   
  if pa_exist==None:
    pa_rt= ResultsTable()   
  else:
    pa_rt = pa_exist.getTextPanel().getOrCreateResultsTable()     
  ParticleAnalyzer.setResultsTable(pa_rt)    
  IJ.run("Set Measurements...", "area mean perimeter shape decimal=3")  
  IJ.run("Analyze Particles...", "size=" + str(cellradius) + "-Infinity circularity=0.1-1.00 add"); #Cell bodies detected 
  pa_rt.show("Results for PA ")
  save_all(srcDir, dstDir, filename, localfile, keepDirectories, imageID)

processFolder()
