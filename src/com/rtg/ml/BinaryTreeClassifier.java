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

import com.rtg.util.MathUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;

/**
 * Implements a binary classification tree.
 */
public final class BinaryTreeClassifier implements PredictClassifier {

  static final int SERIAL_VERSION = 1;

  final BinarySplitter mDirector;
  final PredictClassifier mLeft;
  final PredictClassifier mRight;
  final double mLeftFraction;
  final double mRightFraction;

  /**
   * Constructs a binary tree. When the split criterion cannot decide which subtree to delegate to, a weighted sum is returned.
   * @param director the split criterion
   * @param left the classifier that operates on all instances passing to the left
   * @param right the classifier that operates on all instances passing to the right
   * @param leftFraction the proportion of weight to assign to predictions from the left
   */
  public BinaryTreeClassifier(BinarySplitter director, PredictClassifier left, PredictClassifier right, double leftFraction) {
    if (leftFraction < 0 || leftFraction > 1) {
      throw new IllegalArgumentException("Left fraction must be between 0 and 1");
    }
    if (left == null || right == null) {
      throw new NullPointerException("Left and right trees must not be null");
    }
    mDirector = director;
    mLeft = left;
    mRight = right;
    mLeftFraction = leftFraction;
    mRightFraction = 1.0 - leftFraction;
  }

  /**
   * load a classifier from the stream
   * @param dis the stream to load from
   * @param data set of attributes for encoding/decoding values
   * @throws IOException if an IO error occurs
   */
  public BinaryTreeClassifier(DataInputStream dis, Dataset data) throws IOException {
    final int version = dis.readInt();
    if (version == 1) {
      mDirector = new BinarySplitter(dis, data);
      mLeft = MlPredictLoader.loadPredictClassifier(dis, data);
      mRight = MlPredictLoader.loadPredictClassifier(dis, data);
      mLeftFraction = dis.readDouble();
      mRightFraction = dis.readDouble();
    } else {
      throw new IOException("Unsupported tree version: " + version);
    }
  }

  @Override
  public void save(DataOutputStream dos, Dataset data) throws IOException {
    dos.writeInt(MlPredictLoader.MlPredictType.BINARY_TREE.ordinal()); //tells loader which class to load
    dos.writeInt(SERIAL_VERSION);
    mDirector.save(dos, data);
    mLeft.save(dos, data);
    mRight.save(dos, data);
    dos.writeDouble(mLeftFraction);
    dos.writeDouble(mRightFraction);
  }

  @Override
  public double predict(double[] instance) {
    switch (mDirector.split(instance)) {
      case MISSING:
        return mLeftFraction * mLeft.predict(instance) + mRightFraction * mRight.predict(instance);
      case LEFT:
        return mLeft.predict(instance);
      case RIGHT:
      default:
        return mRight.predict(instance);
    }
  }

  @Override
  public boolean equals(Object obj) {
    if (obj == null || !(obj instanceof BinaryTreeClassifier)) {
      return false;
    }
    final BinaryTreeClassifier obj1 = (BinaryTreeClassifier) obj;
    return mDirector.equals(obj1.mDirector) && mLeft.equals(obj1.mLeft) && mRight.equals(obj1.mRight)
            && MathUtils.approxEquals(mLeftFraction, obj1.mLeftFraction, 10e-10) && MathUtils.approxEquals(mRightFraction, obj1.mRightFraction, 10e-10);
  }

  @Override
  public int hashCode() {
    return Utils.pairHashContinuous(mDirector.hashCode(), mLeft.hashCode(), mRight.hashCode(), Double.valueOf(mLeftFraction).hashCode(), Double.valueOf(mRightFraction).hashCode());
  }

  @Override
  public StringBuilder toString(StringBuilder out, String indent, Dataset data) {
    out.append(indent).append(mDirector.toString(data))
      //.append(" (").append(Utils.realFormat(100.0 * mLeftFraction, 1)).append("%)")
      .append(StringUtils.LS);
    final String newIndent = indent + "  ";
    mLeft.toString(out, newIndent, data);
    mRight.toString(out, newIndent, data);
    return out;
  }
}
