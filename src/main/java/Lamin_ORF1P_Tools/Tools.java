package Lamin_ORF1P_Tools;

import Lamin_ORF1P_Tools.Cellpose.CellposeSegmentImgPlusAdvanced;
import Lamin_ORF1P_Tools.Cellpose.CellposeTaskSettings;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.Duplicator;
import ij.plugin.RGBStackMerge;
import ij.plugin.ZProjector;
import ij.plugin.filter.Analyzer;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.measurements.MeasureCompactness;
import mcib3d.geom2.measurements.MeasureEllipsoid;
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.image3d.ImageHandler;
import org.apache.commons.io.FilenameUtils;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;


/**
 * @author phm
 */
public class Tools {
    public final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    
    public Calibration cal = new Calibration();
    public float pixVol = 0;  
    public String[] chNames = {"Nucleus", "Lamin", "ORF1P"};
    
    // Cellpose
    public int cellPoseDiameter = 100;
    public String cellPoseModel = "cyto2";
    public String cellPoseEnvDirPath = (IJ.isWindows()) ? System.getProperty("user.home")+"\\miniconda3\\envs\\CellPose" : "/opt/miniconda3/envs/cellpose";
    public double minNucVol= 1000;
    public double maxNucVol = 15000; 
    
    private final CLIJ2 clij2 = CLIJ2.getInstance();
 
   
    
    /**
     * Check that needed modules are installed
     */
    public boolean checkInstalledModules() {
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("mcib3d.geom.Object3D");
        } catch (ClassNotFoundException e) {
            IJ.showMessage("Error", "3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        try {
            loader.loadClass("net.haesleinhuepf.clij2.CLIJ2");
        } catch (ClassNotFoundException e) {
            IJ.log("CLIJ not installed, please install from update site");
            return false;
        }
        return true;
    }
    
    
    /**
     * Find images in folder
     */
    public ArrayList findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No image found in " + imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            // Find images with extension
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt) && !f.startsWith("."))
                images.add(imagesFolder + File.separator + f);
        }
        Collections.sort(images);
        return(images);
    }
    
    
    /**
     * Find image calibration
     */
    public Calibration findImageCalib(IMetadata meta) {
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("XY calibration = " + cal.pixelWidth + ", Z calibration = " + cal.pixelDepth);
        return(cal);
    }
    
    
     /**
     * Find channels name
     * @param imageName
     * @return 
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public String[] findChannels (String imageName, IMetadata meta, ImageProcessorReader reader) throws DependencyException, ServiceException, FormatException, IOException {
        ArrayList<String> channels = new ArrayList<>();
        int chs = reader.getSizeC();
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                {
                    if (meta.getChannelID(0, n) == null)
                        channels.add(Integer.toString(n));
                    else 
                        channels.add(meta.getChannelName(0, n).toString());
                }
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null || meta.getChannelName(0, n) == null)
                        channels.add(Integer.toString(n));
                    else 
                        channels.add(meta.getChannelName(0, n).toString());
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels.add(Integer.toString(n));
                    else 
                        channels.add(meta.getChannelFluor(0, n).toString());
                break;
            case "ics" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels.add(Integer.toString(n));
                    else 
                        channels.add(meta.getChannelExcitationWavelength(0, n).value().toString());
                break; 
            case "ics2" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels.add(Integer.toString(n));
                    else 
                        channels.add(meta.getChannelExcitationWavelength(0, n).value().toString());
                break;        
            default :
                for (int n = 0; n < chs; n++)
                    channels.add(Integer.toString(n));

        }
        channels.add("None");
        return(channels.toArray(new String[channels.size()]));         
    }
    
    
    /**
     * Generate dialog box
     */
    public String[] dialog(String[] channels) {
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 80, 0);
        gd.addImage(icon);
          
        
        gd.addMessage("Channels", new Font(Font.MONOSPACED , Font.BOLD, 12), Color.blue);
        int index = 0;
        for (String chName: chNames) {
            gd.addChoice(chName + ": ", channels, channels[index]);
            index++;
        }
        
        gd.addMessage("Nuclei detection", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min nucleus volume (µm3):", minNucVol);
        gd.addNumericField("Max nucleus volume (µm3):", maxNucVol);   
        
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("XY calibration (µm):", cal.pixelWidth);
        gd.addNumericField("Z calibration (µm):", cal.pixelDepth);
        gd.showDialog();
        
        String[] ch = new String[chNames.length];
        for (int i = 0; i < chNames.length; i++)
            ch[i] = gd.getNextChoice();
       
        minNucVol = (float) gd.getNextNumber();
        maxNucVol = (float) gd.getNextNumber();
        
        cal.pixelWidth = cal.pixelHeight = gd.getNextNumber();
        cal.pixelDepth = gd.getNextNumber();
        pixVol = (float) (cal.pixelWidth*cal.pixelHeight*cal.pixelDepth);
        
        if(gd.wasCanceled())
            ch = null;
        
        return(ch);
    }
    
    
    /**
     * Flush and close an image
     */
    public void flush_close(ImagePlus img) {
        img.flush();
        img.close();
    }
    
   
    /**
     * Look for all 3D cells in a Z-stack: 
     * - apply CellPose in 2D slice by slice 
     * - let CellPose reconstruct cells in 3D using the stitch threshold parameters
     */
    public ArrayList<Nucleus> cellposeDetection(ImagePlus img, boolean resize, String cellposeModel, int diameter, double stitchThreshold, boolean zFilter) throws IOException{
        float resizeFactor;
        ImagePlus imgResized;
        if (resize) {
            resizeFactor = 0.5f;
            imgResized = img.resize((int)(img.getWidth()*resizeFactor), (int)(img.getHeight()*resizeFactor), 1, "none");
        } else {
            resizeFactor = 1;
            imgResized = new Duplicator().run(img).duplicate();
        }

        // Define CellPose settings
        CellposeTaskSettings settings = new CellposeTaskSettings(cellposeModel, 1, (int)(diameter*resizeFactor), cellPoseEnvDirPath);
        settings.setStitchThreshold(stitchThreshold);
        settings.useGpu(true);
       
        // Run CellPose
        CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, imgResized);
        ImagePlus imgOut = cellpose.run();
        if(resize) imgOut = imgOut.resize(img.getWidth(), img.getHeight(), "none");
        imgOut.setCalibration(cal);
       
        // Get cells as a population of objects
        ImageHandler imgH = ImageHandler.wrap(imgOut);
        Objects3DIntPopulation pop = new Objects3DIntPopulation(imgH);
        System.out.println(pop.getNbObjects() + " CellPose detections");
       
        // Filter cells by size
        if (zFilter)
            pop = zFilterPop(pop);
        popFilterSize(pop, minNucVol, maxNucVol);
        pop.resetLabels();
        System.out.println(pop.getNbObjects() + " detections remaining after size filtering");
        
        ArrayList<Nucleus> nuclei = new ArrayList<Nucleus>();
        for (Object3DInt obj : pop.getObjects3DInt())
            nuclei.add(new Nucleus(obj));
       
        flush_close(imgOut);
        imgH.closeImagePlus();
        
        return(nuclei);
    } 
    
    
    /*
     * Only keep objects in population that appear in one z-slice only 
     */
    public Objects3DIntPopulation zFilterPop (Objects3DIntPopulation pop) {
        Objects3DIntPopulation popZ = new Objects3DIntPopulation();
        for (Object3DInt obj : pop.getObjects3DInt()) {
            int zmin = obj.getBoundingBox().zmin;
            int zmax = obj.getBoundingBox().zmax;
            if (zmax != zmin)
                popZ.addObject(obj);
        }
        return popZ;
    }
    
    
    /**
     * Filter objects in population by size
     */
    public void popFilterSize(Objects3DIntPopulation pop, double min, double max) {
        pop.getObjects3DInt().removeIf(p -> (new MeasureVolume(p).getVolumeUnit() < min) || (new MeasureVolume(p).getVolumeUnit() > max));
    }
    
    
    /**
     * Compute cells parameters
     */
    public HashMap<String, Double> tagCells(ImagePlus imgORF1P, ImagePlus imgLamin, ImagePlus imgCyto, ArrayList<Nucleus> nuclei) {
        ImageHandler imhORF1P = ImageHandler.wrap(imgORF1P);
        ImageHandler imhLamin = (imgLamin == null) ? null : ImageHandler.wrap(imgLamin);
        // Compute background
        double bgORF1P = findBackground(imgORF1P);
        double bgLamin = (imgLamin == null) ? 0 : findBackground(imgLamin);
        
        // Get cells cytoplasm parameters
        double[] cytoParams = computeCytoParameters(imgCyto, imgORF1P);
        double cytoVol = cytoParams[0];
        double cytoInt = cytoParams[1] - bgORF1P*cytoVol/pixVol;

        // Save cells parameters
        double allNucVol = 0;
        double nucLaminIntSum = 0;
        double nucORF1PIntSum = 0;
        for (Nucleus nucleus: nuclei) {
            Object3DInt nucObj = nucleus.nucleus;

            double nucVol = new MeasureVolume(nucObj).getVolumeUnit();
            allNucVol += nucVol;
            
            double nucComp = new MeasureCompactness(nucObj).getValueMeasurement(MeasureCompactness.COMP_CORRECTED);
            double nucSph = new MeasureCompactness(nucObj).getValueMeasurement(MeasureCompactness.SPHER_CORRECTED);
            double nucElongation = new MeasureEllipsoid(nucObj).getValueMeasurement(MeasureEllipsoid.ELL_ELONGATION);
            double nucFlatness = new MeasureEllipsoid(nucObj).getValueMeasurement(MeasureEllipsoid.ELL_FLATNESS);
            
            double nucORF1PInt = new MeasureIntensity(nucObj, imhORF1P).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
            nucORF1PInt = nucORF1PInt - bgORF1P*nucVol/pixVol;
            nucORF1PIntSum += nucORF1PInt;
            double nucLaminInt = 0;
            if (imhLamin != null) {
                nucLaminInt = new MeasureIntensity(nucObj, imhLamin).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
                nucLaminInt = nucLaminInt - bgLamin*nucVol/pixVol;
                nucLaminIntSum += nucLaminInt;
            }
            
            nucleus.setParams(nucObj.getLabel(), nucVol, nucComp, nucSph, nucElongation, nucFlatness, nucLaminInt, nucORF1PInt);
        }
        // Save global parameters
        HashMap<String, Double> globalParams = new HashMap<>();
        globalParams.put("nucNb", (double) nuclei.size());
        globalParams.put("nucVol", allNucVol);
        globalParams.put("bgLamin", bgLamin);
        globalParams.put("nucLaminIntSum", nucLaminIntSum);
        globalParams.put("bgORF1P", bgORF1P);
        globalParams.put("nucORF1PIntSum", nucORF1PIntSum);
        globalParams.put("cytoORF1PVol", cytoVol);
        globalParams.put("cytoORF1PInt", cytoInt);
        
        imhORF1P.closeImagePlus();
        if (imhLamin != null)
            imhLamin.closeImagePlus();

        return(globalParams);
    }
    
    
    /**
     * Find background image intensity:
     * Z projection over min intensity + read median intensity
     */
    public double findBackground(ImagePlus img) {
      ImagePlus imgProj = doZProjection(img, ZProjector.MIN_METHOD);
      ImageProcessor imp = imgProj.getProcessor();
      double bg = imp.getStatistics().median;
      flush_close(imgProj);
      return(bg);
    }
    
    
    /**
     * Do Z projection
     */
    public ImagePlus doZProjection(ImagePlus img, int param) {
        ZProjector zproject = new ZProjector();
        zproject.setMethod(param);
        zproject.setStartSlice(1);
        zproject.setStopSlice(img.getNSlices());
        zproject.setImage(img);
        zproject.doProjection();
       return(zproject.getProjection());
    }
    
    
    /**
     * Detect cytoplasm and compute its volume
     */    
    public ImagePlus findCellsCyto(ImagePlus img, ArrayList<Nucleus> nuclei) {
        ImagePlus imgIn = new Duplicator().run(img).duplicate();
        ImagePlus imgTh = gaussian_filter(imgIn, 4, 4);
        imgTh.setSlice(imgTh.getNSlices()/2);
        IJ.setAutoThreshold(imgTh, "Huang dark");
        Prefs.blackBackground = false;
        IJ.run(imgTh, "Convert to Mask", "method=Huang background=Dark");
        IJ.run(imgTh, "Median...", "radius=8 stack");

        ImageHandler imh = ImageHandler.wrap(imgTh);
        for (Nucleus nucleus : nuclei)
            nucleus.nucleus.drawObject(imh, 0);
        imgTh = imh.getImagePlus();

        flush_close(imgIn);
        return(imgTh);
    }
    
    
    /**
     * Gaussian filter using CLIJ2
     */ 
    public ImagePlus gaussian_filter(ImagePlus img, double sizeXY, double sizeZ) {
       ClearCLBuffer imgCL = clij2.push(img);
       ClearCLBuffer imgCLG = clij2.create(imgCL);
       clij2.gaussianBlur3D(imgCL, imgCLG, sizeXY, sizeXY, sizeZ);
       ImagePlus imgG = clij2.pull(imgCLG);
       imgG.setCalibration(cal);
       clij2.release(imgCL);
       clij2.release(imgCLG);
       return(imgG);
    } 
    
    
    public double[] computeCytoParameters(ImagePlus mask, ImagePlus img) {
        ResultsTable rt = new ResultsTable();
        Analyzer analyzer = new Analyzer(img, Analyzer.AREA+Analyzer.INTEGRATED_DENSITY, rt);
        double area = 0;
        double intSum = 0;
        for (int n = 1; n <= mask.getNSlices(); n++) {
            mask.setSlice(n);
            IJ.setAutoThreshold(mask, "Default");
            IJ.run(mask, "Create Selection", "");
            Roi roi = mask.getRoi();
            
            img.setSlice(n);
            img.setRoi(roi);
            rt.reset();
            analyzer.measure();
            area += rt.getValue("Area", 0); 
            intSum += rt.getValue("RawIntDen",0);
        }
        double[] volInt = {area*cal.pixelDepth, intSum};
        return(volInt);
    }
    
    
    /**
     * Draw results in images
     */
    public void drawResults(ArrayList<Nucleus> nuclei, ImagePlus imgCyto, ImagePlus img, String imgName, String outDir) {
        ImageHandler imgNuc = ImageHandler.wrap(new Duplicator().run(img)).createSameDimensions();
        for (Nucleus nucleus: nuclei)
            nucleus.nucleus.drawObject(imgNuc);
        new ImageConverter(imgCyto).convertToGray16();
        
        ImagePlus[] imgColors = {null, imgCyto, imgNuc.getImagePlus(), img};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, true);
        imgObjects.setCalibration(cal);
        FileSaver ImgObjectsFile1 = new FileSaver(imgObjects);
        ImgObjectsFile1.saveAsTiff(outDir + imgName + ".tif");
        
        flush_close(imgObjects);
        imgNuc.closeImagePlus();
    }
    
}