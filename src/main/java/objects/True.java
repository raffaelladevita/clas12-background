package objects;

import org.jlab.detector.base.DetectorType;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;
import org.jlab.io.base.DataBank;

/**
 *
 * @author devita
 */
public class True implements Comparable {
 
    private final int hitn;
    private final DetectorType type;
    
    private double time;
    private double edep;
    private double energy;
    private Vector3D momentum;
    private Point3D  position;
    private Point3D  vertex;
    private int tid;
    private int pid;
    private int mpid;
    
    
    public True(int hitn, DetectorType type) {
        this.hitn = hitn;
        this.type = type;
    }

    public DetectorType getType() {
        return type;
    }

    public int getHitn() {        
        return hitn;
    }    
    
    public double getKinEnergy() {        
        return energy-Math.sqrt(Math.max(0, energy*energy-momentum.mag2()));
    }

    public double getEnergy() {        
        return energy;
    }

    public void setEnergy(double energy) {
        this.energy = energy;
    }

    public double getTime() {
        return time;
    }

    public void setTime(double time) {
        this.time = time;
    }

    public double getEdep() {
        return edep;
    }

    public void setEdep(double edep) {
        this.edep = edep;
    }

    public Vector3D getMomentum() {
        return momentum;
    }

    public void setMomentum(Vector3D momentum) {
        this.momentum = momentum;
    }

    public Point3D getPosition() {
        return position;
    }

    public void setPosition(Point3D position) {
        this.position = position;
    }

    public Point3D getVertex() {
        return vertex;
    }

    public void setVertex(Point3D vertex) {
        this.vertex = vertex;
    }

    public int getTid() {
        return tid;
    }

    public void setTID(int tid) {
        this.tid = tid;
    }

    public int getPid() {
        return pid;
    }

    public void setPID(int pid) {
        this.pid = pid;
    }

    public int getMPID() {
        return mpid;
    }

    public void setMPID(int mpid) {
        this.mpid = mpid;
    }
    
    public String getName() {
        return this.type.getName();
    }
    
    public static True readTruth(DataBank mc, int row) {
        True hit = new True(mc.getInt("hitn", row),
                            DetectorType.getType(mc.getInt("detector", row)));
        hit.setEnergy(mc.getFloat("trackE", row));
        hit.setEdep(mc.getFloat("totEdep", row));
        hit.setTime(mc.getFloat("avgT", row));
        hit.setMomentum(new Vector3D(mc.getFloat("px", row),
                                     mc.getFloat("py", row),
                                     mc.getFloat("pz", row)));
        hit.setPosition(new Point3D(mc.getFloat("avgX", row),
                                    mc.getFloat("avgY", row),
                                    mc.getFloat("avgZ", row)));
        hit.setVertex(new Point3D(mc.getFloat("vx", row),
                                    mc.getFloat("vy", row),
                                    mc.getFloat("vz", row)));
        hit.setTID(mc.getInt("tid", row));
        hit.setPID(mc.getInt("pid", row));
        hit.setMPID(mc.getInt("mpid", row));
        return hit;
    }

    @Override
    public int compareTo(Object o) {
        True ot = (True) o;
        if(this.getTime()<ot.getTime()) return -1;
        else if(this.getTime()==ot.getTime()) return 0;
        else return 1;
    }
    
    @Override
    public String toString() {
        String s = String.format("True: detector=%s id=%d pid=%d time=%.3f ns E=%.3f MeV p=%.3f MeV, Ekin=%.3f MeV", 
                                  this.type.getName(), this.hitn, this.pid, this.time, this.energy, this.momentum.mag(), this.getKinEnergy());
        return s;
    }
    
}
