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

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.Utils;

/**
 * This class evaluates an instance with respect to a binary decision. For numeric attributes, values that
 * are "less than or equal to" the given split value proceed left. For nominal attributes, values that are equal to
 * the split value proceed left.
 */
@TestClass(value = {"com.rtg.ml.BinaryTreeClassifierTest", "com.rtg.ml.BinarySplitterTest"})
public final class BinarySplitter {

  static final int SERIAL_VERSION = 2;

  /**
   * Encapsulates the directions that can be taken after a split
   */
  public enum Direction {
    /** Go left */
    LEFT,
    /** Go right */
    RIGHT,
    /** Cannot determine direction */
    MISSING
  }

  private final String mName;
  private final int mAttributeIndex;
  private final MlDataType mSplitValueDataType;
  private final double mSplitValue;
  private final boolean mNumeric;
  private final boolean mSplitMissing;

  /**
   * Constructor
   * @param name the name of the attribute
   * @param attribute the index of the attribute to consider
   * @param splitValue the split value.
   * @param splitValueDataType data type of {@code splitValue}
   */
  public BinarySplitter(String name, int attribute, double splitValue, MlDataType splitValueDataType) {
    mName = name;
    mAttributeIndex = attribute;
    mSplitValueDataType = splitValueDataType;
    mSplitValue = splitValue;
    mSplitMissing = Attribute.isMissingValue(splitValue);
    mNumeric = splitValueDataType.isNumeric();
  }

  /**
   * @param dis stream to load splitter from
   * @param data set of attributes for encoding/decoding values
   * @throws IOException if an IO error occurs, or a newer version is attempted to be read
   */
  public BinarySplitter(DataInputStream dis, Dataset data) throws IOException {
    final int version = dis.readInt();

    if (version == 1) {
      dis.readUTF(); // Ignored, instead get name from dataset and attribute index
      mAttributeIndex = dis.readInt();
      mName = data.getAttributes()[mAttributeIndex].getName();
      final int type = dis.readInt();
      if (type >= MlDataType.values().length || type < 0) {
        throw new IOException("Learning attribute out of range, model may be corrupt or created with a later version of RTG.");
      }
      mSplitValueDataType = MlDataType.values()[type];
      mSplitValue = data.getAttributes()[mAttributeIndex].encodeValue(mSplitValueDataType.load(dis));
      mSplitMissing = false;
      mNumeric = dis.readBoolean();

    } else if (version == 2) {
      mAttributeIndex = dis.readInt();
      mName = data.getAttributes()[mAttributeIndex].getName();
      final int type = dis.readInt();
      if (type >= MlDataType.values().length || type < 0) {
        throw new IOException("Learning attribute out of range, model may be corrupt or created with a later version of RTG.");
      }
      mSplitValueDataType = MlDataType.values()[type];
      mSplitMissing = dis.readBoolean();
      if (mSplitMissing) {
        mSplitValue = Double.NaN;
      } else {
        mSplitValue = data.getAttributes()[mAttributeIndex].encodeValue(mSplitValueDataType.load(dis));
      }
      mNumeric = mSplitValueDataType.isNumeric();

    } else {
      throw new IOException("Unsupported binary split version, model may be corrupt or created with a later version of RTG.");
    }
  }

  /**
   * save the splitter
   * @param dos the stream to save to
   * @param data set of attributes for encoding/decoding values
   * @throws IOException if an IO error occurs
   */
  public void save(DataOutputStream dos, Dataset data) throws IOException {
    save(dos, data, SERIAL_VERSION);
  }

  void save(DataOutputStream dos, Dataset data, int version) throws IOException {
    switch (version) {
      case 1:
        saveV1(dos, data);
        break;
      case 2:
        saveV2(dos, data);
        break;
      default:
        throw new IOException("Can not save as version " + version);
    }
  }

  void saveV1(DataOutputStream dos, Dataset data) throws IOException {
    if (mSplitMissing) {
      throw new IOException("Cannot save split-on-missing in version 1");
    }
    final Attribute attribute = data.getAttributes()[mAttributeIndex];
    dos.writeInt(1);
    dos.writeUTF(attribute.getName());
    dos.writeInt(mAttributeIndex);
    dos.writeInt(mSplitValueDataType.ordinal());
    mSplitValueDataType.save(attribute.decodeValue(mSplitValue), dos);
    dos.writeBoolean(mNumeric);
  }

  void saveV2(DataOutputStream dos, Dataset data) throws IOException {
    dos.writeInt(2);
    dos.writeInt(mAttributeIndex);
    dos.writeInt(mSplitValueDataType.ordinal());
    dos.writeBoolean(mSplitMissing);
    if (!mSplitMissing) {
      final Attribute attribute = data.getAttributes()[mAttributeIndex];
      mSplitValueDataType.save(attribute.decodeValue(mSplitValue), dos);
    }
  }

  /**
   * Evaluate which path of a binary tree this instance should take
   * @param instance the instance
   * @return LEFT if the instance should go left, RIGHT if it should go right, or MISSING if the direction could
   * not be determined
   */
  public Direction split(double[] instance) {
    return split(instance[mAttributeIndex]);
  }

  Direction split(double splitValue) {
    if (mSplitMissing) {
      return Attribute.isMissingValue(splitValue) ? Direction.LEFT : Direction.RIGHT;
    } else if (Attribute.isMissingValue(splitValue)) {
      return Direction.MISSING;
    } else if (mNumeric) { // Numeric
      try {
        @SuppressWarnings(value = {"unchecked", "rawtypes"})
        final Direction retVal = ((Comparable) mSplitValue).compareTo(splitValue) >= 0 ? Direction.LEFT : Direction.RIGHT;
        return retVal;
      } catch (ClassCastException e) {
        throw new RuntimeException("Illegal value type for attribute " + mName + ".", e);
      }
    } else { // Other nominal attributes take a binary split
      //in this case the double is effectively an enumeration, so should only take integer values
      return splitValue == mSplitValue ? Direction.LEFT : Direction.RIGHT;
    }
  }

  @Override
  public int hashCode() {
    return Utils.pairHashContinuous(mName.hashCode(), mAttributeIndex, mSplitMissing ? 42 : Double.valueOf(mSplitValue).hashCode(), mSplitValueDataType.hashCode(), Boolean.valueOf(mNumeric).hashCode());
  }

  @Override
  public boolean equals(Object obj) {
    if (obj == null || !(obj instanceof  BinarySplitter)) {
      return false;
    }
    final BinarySplitter obj1 = (BinarySplitter) obj;
    return mName.equals(obj1.mName)
      && mAttributeIndex == obj1.mAttributeIndex
      && mSplitValueDataType == obj1.mSplitValueDataType
      && mNumeric == obj1.mNumeric
      && (mSplitValue == obj1.mSplitValue || (mSplitMissing && obj1.mSplitMissing));
  }

  /**
   * @param data set of attributes for encoding/decoding values
   * @return string representation of split point
   */
  public String toString(Dataset data) {
    final StringBuilder out = new StringBuilder();
    out.append("split: ").append(data.getAttributes()[mAttributeIndex].getName());
    if (mSplitMissing) {
      out.append(" == missing");
    } else if (mNumeric) {
      out.append(" <= ").append(data.getAttributes()[mAttributeIndex].decodeValue(mSplitValue));
    } else {
      out.append(" == ").append(data.getAttributes()[mAttributeIndex].decodeValue(mSplitValue));
    }
    return out.toString();
  }
}
