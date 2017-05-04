package mdsim;
import java.io.*;
import java.util.*;


enum Arithmetic {
  ADD, SUB, MUL, DIV
}

class Coordinate implements Serializable,Comparable {
  static final int DIMENSION = 3;
  public double[] component;

  public Coordinate() {
    component = new double[DIMENSION];
  }

  public Coordinate(Coordinate val) {
    component = Arrays.copyOf(val.component, val.component.length);
  }

  public Coordinate(double v) {
    this();
    Arrays.fill(component, v);
  }

  public Coordinate(double[] val) {
    component = val.clone();
  }
  
  public void set(double[] val) {
    component = val.clone();
  }

  public Coordinate get() {
    return new Coordinate(component);
  }

  public Coordinate operation(Object o, Arithmetic op) {
    double v2[];
    if (o instanceof Double) {
      double v = ((Double)o);
      v2 = new double[component.length];
      Arrays.fill(v2, v);
    }
    else {
      v2 = ((Coordinate)o).component;
    }

    for (int i=0; i<component.length; ++i) {
      switch (op) {
        case ADD:
          component[i] += v2[i];
          break;
        case SUB:
          component[i] -= v2[i];
          break;
        case MUL:
          component[i] *= v2[i];
          break;
        case DIV:
          component[i] /= v2[i];
          break;
      }
    }
    return new Coordinate(component);
  }

  @Override
  public String toString() {
    return String.format("%f %f %f", component[0], component[1], component[2]);
  }

  @Override                                                                           
  public int compareTo(Object o) {                                                    
    Coordinate s = (Coordinate)o;
    for (int i=0; i<this.component.length; ++i) {
      if (this.component[i] != s.component[i]) {
        return Double.compare(this.component[i], s.component[i]);
      }
    }
    return 0;
  }
}


class Particle extends Coordinate {
  static final double LJEplison=1.0, LJSigma=1.0, Mass=1.0;
  static final double LJSigma6 = LJSigma*LJSigma*LJSigma*LJSigma*LJSigma*LJSigma;
  static final double BoxSize=15.0;
  public double old[];

  public Particle(double[] val) {
    super(val);
    old = Arrays.copyOf(component, component.length);
  }

  public Particle(Particle val) {
    super(val);
    old = Arrays.copyOf(val.old, val.old.length);
  }

  public void set(Coordinate val) {
    set(val.component);
  }

  @Override
  public void set(double[] val) {
    old = component;
    super.set(val);
    for (int i=0; i<component.length; ++i) {
      if (component[i] >  BoxSize) {
        component[i] -= 2*BoxSize;
        old[i]       -= 2*BoxSize;
      }
      else if (component[i] < -BoxSize) {
        component[i] += 2*BoxSize;
        old[i]       += 2*BoxSize;
      }
    }
  }
  
  public Coordinate get_old() {
    return new Coordinate(old);
  }

  public double potentialTo(Particle s) {
    double dist2 = 0;
    for (int i=0; i<this.component.length; ++i) {
      double d = this.component[i] - s.component[i];
      if      (d >  BoxSize) d -= 2*BoxSize;
      else if (d < -BoxSize) d += 2*BoxSize;
      dist2 += d*d;
    }
    double dist6 = dist2*dist2*dist2;
    double sigdist6 = LJSigma6 / dist6;
    return 4 * LJEplison * (sigdist6*sigdist6 - sigdist6);
  }
}

class WorldRange implements Serializable {
  public Particle range[];
  public int range_begin;
  public WorldRange(Particle[] world, int range_begin, int range_len) {
    this.range = Arrays.copyOfRange(world, range_begin, range_begin+range_len);
    this.range_begin = range_begin;
  }
}


