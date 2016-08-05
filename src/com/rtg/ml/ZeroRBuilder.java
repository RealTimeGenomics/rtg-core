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
    static final int SERIAL_VERSION = 1;
    private final double mProb;
    private final String mDesc;
    final int mCurrentVersion;

    /**
     * @param prob weight to give to positive
     */
    public ZeroRClassifier(double prob) {
      mProb = prob;
      mDesc = "" + prob;
      mCurrentVersion = SERIAL_VERSION;
    }

    /**
     * @param pos number (or total weight) of positive examples
     * @param neg number (or total weight) of negative examples
     */
    public ZeroRClassifier(double pos, double neg) {
      if (MathUtils.approxEquals(pos + neg, 0, MINIMUM_WEIGHT)) {
        throw new IllegalArgumentException("0R undefined on no examples");
      }
      mProb = pos / (pos + neg);
      mDesc = "" + mProb + " (" + pos + "/" + (pos + neg) + ")";
      mCurrentVersion = SERIAL_VERSION;
    }

    /**
     * @param dis stream to load from
     * @throws IOException if an IO error occurs
     */
    public ZeroRClassifier(DataInputStream dis) throws IOException {
      mCurrentVersion = dis.readInt();
      if (mCurrentVersion == 1) {
        mProb = dis.readDouble();
        mDesc = dis.readUTF();
      } else {
        throw new IOException("Unsupported 0R version: " + mCurrentVersion);
      }
    }

    @Override
    public void save(DataOutputStream dos, Dataset data) throws IOException {
      dos.writeInt(MlPredictLoader.MlPredictType.ZERO_R.ordinal()); //tell loader which class to load from
      dos.writeInt(SERIAL_VERSION);
      dos.writeDouble(mProb);
      dos.writeUTF(mDesc);
    }

    @Override
    public boolean equals(Object o) {
      if (o == null || !(o instanceof ZeroRClassifier)) {
        return false;
      }
      final ZeroRClassifier other = (ZeroRClassifier) o;
      return MathUtils.approxEquals(mProb, other.mProb, 10e-10) && mDesc.equals(other.mDesc);
    }

    @Override
    public int hashCode() {
      return Utils.pairHash(Double.valueOf(mProb).hashCode(), mDesc.hashCode());
    }

    @Override
    public double predict(double[] instance) {
      return mProb;
    }

    @Override
    public StringBuilder toString(StringBuilder out, String indent, Dataset data) {
      return out.append(indent).append("0R: ").append(mDesc).append(StringUtils.LS);
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
