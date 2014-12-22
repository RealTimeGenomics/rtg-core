/*
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.variant.sv.discord;

import com.rtg.util.machine.MachineOrientation;

import net.sf.samtools.SAMRecord;


/**
 * Describes the orientation on the genome of the constraints of a breakpoint.
 * U(up) means the breakpoint is forward of a specified position on the genome
 * D(down) means the breakpoint is back from a specified position on the genome.
 * The terms forward and reverse are avoided because they are used for read end orientations.
 * The first constraint is taken from the first of the reads (even if it is not the actual
 * read seen at this point).
 * The first letter refers to the read arm being referred to and the second letter its mate. This means
 * that a read arm described as <code>UD</code> will have a mate described as <code>DU</code>.
 */
public enum Orientation {
  /** Up Up*/
  UU(+1, +1) {
    @Override
    public Orientation flip() {
      return UU;
    }
    @Override
    public int x(int x) {
      return x;
    }
    @Override
    public int y(int y) {
      return y;
    }
    @Override
    public double x(double x) {
      return x;
    }
    @Override
    public double y(double y) {
      return y;
    }
  },
  /** Up Down */
  UD(+1, -1) {
    @Override
    public Orientation flip() {
      return DU;
    }
    @Override
    public int x(int x) {
      return x;
    }
    @Override
    public int y(int y) {
      return -y;
    }
    @Override
    public double x(double x) {
      return x;
    }
    @Override
    public double y(double y) {
      return -y;
    }
  },
  /** Down Up */
  DU(-1, +1) {
    @Override
    public Orientation flip() {
      return UD;
    }
    @Override
    public int x(int x) {
      return -x;
    }
    @Override
    public int y(int y) {
      return y;
    }
    @Override
    public double x(double x) {
      return -x;
    }
    @Override
    public double y(double y) {
      return y;
    }
  },
  /** Down Down */
  DD(-1, -1)  {
    @Override
    public Orientation flip() {
      return DD;
    }
    @Override
    public int x(int x) {
      return -x;
    }
    @Override
    public int y(int y) {
      return -y;
    }
    @Override
    public double x(double x) {
      return -x;
    }
    @Override
    public double y(double y) {
      return -y;
    }
  };

  private final int mX;
  private final int mY;
  /**
   * @param x direction of the first read
   * @param y direction of the second read
   */
  private Orientation(int x, int y) {
    mX = x;
    mY = y;
  }
  /**
   * Get x.
   * @return Returns the x where +1 means up and -1 means down
   */
  public int getX() {
    return mX;
  }
  /**
   * Get y.
   * @return Returns the y where +1 means up and -1 means down
   */
  public int getY() {
    return mY;
  }

  /**
   * @return the orientation with x and y axes interchanged.
   */
  public abstract Orientation flip();

  /**
   * Transform along the x-axis.
   * @param x value to be transformed.
   * @return the transformed value.
   */
  public abstract int x(final int x);

  /**
   * Transform along the ys-axis.
   * @param y value to be transformed.
   * @return the transformed value.
   */
  public abstract int y(final int y);

  /**
   * Transform along the x-axis.
   * @param x value to be transformed.
   * @return the transformed value.
   */
  public abstract double x(final double x);

  /**
   * Transform along the ys-axis.
   * @param y value to be transformed.
   * @return the transformed value.
   */
  public abstract double y(final double y);

  /**
   * Get orientation given differences between x and y co-ordinates.
   * @param xDir usually computed as <code>z-x</code>
   * @param yDir usually computed as <code>w-y</code>
   * @return the orientation.
   */
  public static Orientation orientation(final int xDir, final int yDir) {
    //In real data shouldnt get 0 for either but treat as positive.
    if (xDir < 0) {
      if (yDir < 0) {
        return DD;
      } else {
        return DU;
      }
    } else {
      if (yDir < 0) {
        return UD;
      } else {
        return UU;
      }
    }
  }

  /**
   * @param rec the sam record to determine the orientation of
   * @param mo the expected orientation of read arms for the machine type
   * @return the orientation of the sam record
   */
  public static Orientation orientation(SAMRecord rec, MachineOrientation mo) {
    return mo.orientation(rec);
  }

}
