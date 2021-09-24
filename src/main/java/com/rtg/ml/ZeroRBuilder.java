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
package com.rtg.ml;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.Properties;

import com.rtg.util.MathUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;

/**
 * Build a simple classifier that always predicts the most frequent class.
 *
 */
public class ZeroRBuilder implements BuildClassifier {

  private PredictClassifier mClassifier = null;

  /**
   * a simple classifier that always predicts the most frequent class
   */
  public static final class ZeroRClassifier implements PredictClassifier {
    static final double MINIMUM_WEIGHT = 1e-10;
    static final int SERIAL_VERSION = 2;

    private final double mProb;
    private final double mTotal;
    private final String mDesc;

    /**
     * @param pos number (or total weight) of positive examples
     * @param neg number (or total weight) of negative examples
     */
    public ZeroRClassifier(double pos, double neg) {
      if (MathUtils.approxEquals(pos + neg, 0, MINIMUM_WEIGHT)) {
        throw new IllegalArgumentException("0R undefined on no examples");
      }
      mProb = pos / (pos + neg);
      mDesc = null;
      mTotal = pos + neg;
    }

    /**
     * @param dis stream to load from
     * @throws IOException if an IO error occurs
     */
    public ZeroRClassifier(DataInputStream dis) throws IOException {
      final int version = dis.readInt();
      if (version == 1) {
        mProb = dis.readDouble();
        mDesc = dis.readUTF();
        mTotal = Double.NaN;
      } else if (version == 2) {
        mProb = dis.readDouble();
        mDesc = null;
        mTotal = dis.readDouble();
      } else {
        throw new IOException("Unsupported 0R version: " + version);
      }
    }

    @Override
    public void save(DataOutputStream dos, Dataset data) throws IOException {
      save(dos, SERIAL_VERSION);
    }

    void save(DataOutputStream dos, int version) throws IOException {
      switch (version) {
        case 1:
          saveV1(dos);
          break;
        case 2:
          saveV2(dos);
          break;
        default:
          throw new IOException("Can not save as version " + version);
      }
    }

    void saveV1(DataOutputStream dos) throws IOException {
      dos.writeInt(MlPredictLoader.MlPredictType.ZERO_R.ordinal()); //tell loader which class to load from
      dos.writeInt(SERIAL_VERSION);
      dos.writeDouble(mProb);
      dos.writeUTF(mDesc != null ? mDesc : desc());
    }

    void saveV2(DataOutputStream dos) throws IOException {
      dos.writeInt(MlPredictLoader.MlPredictType.ZERO_R.ordinal()); //tell loader which class to load from
      dos.writeInt(SERIAL_VERSION);
      dos.writeDouble(mProb);
      dos.writeDouble(mTotal);
    }

    @Override
    public boolean equals(Object o) {
      if (o == null || !(o instanceof ZeroRClassifier)) {
        return false;
      }
      final ZeroRClassifier other = (ZeroRClassifier) o;
      return MathUtils.approxEquals(mProb, other.mProb, 10e-10);
    }

    @Override
    public int hashCode() {
      return Double.valueOf(mProb).hashCode();
    }

    @Override
    public double predict(double[] instance) {
      return mProb;
    }

    @Override
    public StringBuilder toString(StringBuilder out, String indent, Dataset data) {
      out.append(indent).append("0R: ");
      if (mDesc != null) {
        out.append(mDesc);
      } else {
        out.append(desc());
      }
      out.append(StringUtils.LS);
      return out;
    }
    @Override
    public String toString() {
      return toString(new StringBuilder(), "", null).toString();
    }

    private String desc() {
      return Utils.realFormat(mProb, 4)
        + (Double.isNaN(mTotal) ? "" : String.format(" (%d/%d)", (int) (mProb * mTotal), (int) mTotal));
    }
  }

  @Override
  public void build(Dataset dataset) {
    mClassifier = new ZeroRClassifier(dataset.totalPositiveWeight(), dataset.totalNegativeWeight());
  }

  @Override
  public PredictClassifier getClassifier() {
    return mClassifier;
  }

  @Override
  public void setProperties(Properties props) {
  }

}
