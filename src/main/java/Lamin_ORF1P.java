import Lamin_ORF1P_Tools.Tools;
import Lamin_ORF1P_Tools.Nucleus;
import ij.IJ;
import ij.ImagePlus;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.HashMap;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.util.ImageProcessorReader;
import ij.plugin.PlugIn;
import java.io.FileWriter;
import java.util.ArrayList;
import loci.plugins.in.ImporterOptions;
import org.apache.commons.io.FilenameUtils;
import org.scijava.util.ArrayUtils;


/**
 * Find ORF1P cells and their corresponding DAPI nuclei
 * Measure nucleus and cytoplasm intensity in ORF1P channel (488)
 *
 * @author phm
 */
public class Lamin_ORF1P implements PlugIn {
    
    Tools tools = new Tools();
    private String imageDir = "";
    public  String outDirResults = "";
    public  String rootName = "";
    public BufferedWriter results, globalResults;
    
    public void run(String arg) {
        try {
            if ((!tools.checkInstalledModules())) {
                return;
            } 
      
            imageDir = IJ.getDirectory("Choose directory containing image files...");
            if (imageDir == null) {
                return;
            }
            
            // Find images with nd extension
            ArrayList<String> imageFile = tools.findImages(imageDir, "nd");
            if (imageFile == null) {
                IJ.showMessage("Error", "No images found with nd extension");
                return;
            }
            
            // Create output folder
            outDirResults = imageDir + File.separator + "Results" + File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            
            // Write header for nuclei parameters file
            String header = "Image name\tNucleus ID\tNucleus volume (µm3)\tNucleus compactness\tNucleus sphericity\tNucleus elongation\t"
                    + "Nucleus flatness\tNucleus bg corr. total intensity in lamin channel\tNucleus bg corr. total intensity in ORF1P channel\n";
            FileWriter fwResults = new FileWriter(outDirResults + "detailed_results.xls", false);
            results = new BufferedWriter(fwResults);
            results.write(header);
            results.flush();
            
            // Write header for global nuclei parameters file
            header = "Image name\tNb nuclei\tNuclei total volume (µm3)\tLamin background\tNuclei bg corr. total intensity in lamin channel\t"
                    + "ORF1P background\tNuclei bg corr. total intensity in ORF1P channel\t"
                    + "ORF1P cytoplasm total volume (µm3)\tORF1P cytoplasm bg corr. total intensity\n";
            FileWriter fwGlobalResults = new FileWriter(outDirResults + "global_results.xls", false);
            globalResults = new BufferedWriter(fwGlobalResults);
            globalResults.write(header);
            globalResults.flush();
            
            // Create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imageFile.get(0));
            
            // Find image calibration
            tools.cal = tools.findImageCalib(meta);
            
            // Find channel names
            String[] channels = tools.findChannels(imageFile.get(0), meta, reader);

            // Channels dialog
            String[] chs = tools.dialog(channels);
            if (chs == null) {
                IJ.showStatus("Plugin canceled");
                return;
            }
            
            for (String f : imageFile) {
                rootName = FilenameUtils.getBaseName(f);
                System.out.println("--- ANALYZING IMAGE " + rootName + " ------");
                reader.setId(f);
                
                ImporterOptions options = new ImporterOptions();
                options.setId(f);
                options.setSplitChannels(true);
                options.setQuiet(true);
                options.setCrop(true);
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);

                // Open DAPI channel
                System.out.println("- Analyzing " + chs[0] + " nuclei channel -");
                int indexCh = ArrayUtils.indexOf(channels, chs[0]);
                ImagePlus imgNucleus = BF.openImagePlus(options)[indexCh];
                
                // Find DAPI nuclei
                System.out.println("Finding nuclei....");
                ArrayList<Nucleus> nuclei = tools.cellposeDetection(imgNucleus, true, "cyto2", 100, 0.5, true);
                System.out.println(nuclei.size() + " " + chs[0] + " nuclei found");
                tools.flush_close(imgNucleus);
                
                // Open ORF1P channel
                System.out.println("- Analyzing " + chs[2] + " ORF1P channel -");
                indexCh = ArrayUtils.indexOf(channels, chs[2]);
                ImagePlus imgORF1P = BF.openImagePlus(options)[indexCh];
                
                // Threshold cytoplasm
                System.out.println("Finding cytoplasm....");
                ImagePlus imgCyto = tools.findCellsCyto(imgORF1P, nuclei);
                
                // Open Lamin channel (if provided)
                ImagePlus imgLamin = null;
                if (!chs[1].equals("None")) {
                    System.out.println("- Analyzing " + chs[1] + " lamin channel -");
                    indexCh = ArrayUtils.indexOf(channels, chs[1]);
                    imgLamin = BF.openImagePlus(options)[indexCh];
                }
               
                // Tag nuclei with parameters
                System.out.println("- Measuring cells parameters -");
                HashMap<String, Double> globalParams = tools.tagCells(imgORF1P, imgLamin, imgCyto, nuclei);
                
                // Draw results image
                System.out.println("- Saving results -");
                tools.drawResults(nuclei, imgCyto, imgORF1P, rootName, outDirResults);
                             
                // Write nuclei parameters results
                for (Nucleus nucleus: nuclei) {
                    HashMap<String, Double> params = nucleus.params;
                    results.write(rootName+"\t"+params.get("nucIndex")+"\t"+params.get("nucVol")+"\t"+params.get("nucComp")+
                            "\t"+params.get("nucSph")+"\t"+params.get("nucEllElong")+"\t"+params.get("nucEllFlat")+
                            "\t"+params.get("nucLaminInt")+"\t"+params.get("nucORF1PInt")+"\n");
                    results.flush();
                }
                
                // Write global parameters results
                globalResults.write(rootName+"\t"+globalParams.get("nucNb").intValue()+"\t"+globalParams.get("nucVol")+"\t"+globalParams.get("bgLamin")+
                        "\t"+globalParams.get("nucLaminIntSum")+"\t"+globalParams.get("bgORF1P")+"\t"+globalParams.get("nucORF1PIntSum")+"\t"+globalParams.get("cytoORF1PVol")+"\t"+globalParams.get("cytoORF1PInt")+"\n");
                globalResults.flush();
                
                if (imgLamin != null)
                    tools.flush_close(imgLamin);
                tools.flush_close(imgORF1P);
                tools.flush_close(imgCyto);
            }            
            results.close();
            globalResults.close();
        } catch (IOException | DependencyException | ServiceException | FormatException  ex) {
            Logger.getLogger(Lamin_ORF1P.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        System.out.println("--- All done! ---");
    }
}

           