package Lamin_Tools;

import java.util.HashMap;
import mcib3d.geom2.Object3DInt;

/**
 * @author hm
 */
public class Nucleus {
    
    public Object3DInt nucleus;
    public HashMap<String, Double> params;
    
    public Nucleus(Object3DInt nucleus) {
        this.nucleus = nucleus;
        this.params = new HashMap<>();
    }
    
    
    public void setParams(double nucLabel, double nucVol, double nucComp, double nucSph, double nucEllElong, double nucEllFlat, double nucLaminInt,
            double nucIntORF1P) {
        params.put("nucIndex", nucLabel);
        params.put("nucVol", nucVol);
        params.put("nucComp", nucComp);
        params.put("nucSph", nucSph);
        params.put("nucEllElong", nucEllElong);
        params.put("nucEllFlat", nucEllFlat);
        params.put("nucLaminInt", nucLaminInt);
        params.put("nucORF1PInt", nucIntORF1P);
    }
}
