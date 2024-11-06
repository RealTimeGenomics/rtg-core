/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package com.rtg.variant.avr;

import static com.rtg.util.StringUtils.LS;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.SortedSet;
import java.util.TreeSet;

import com.rtg.ml.Attribute;
import com.rtg.ml.Dataset;
import com.rtg.ml.MlDataType;
import com.rtg.util.TextTable;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfHeader;

/**
 * Provides conversion from VCF records into instances for the machine learning.
 * Canonicalization of annotations/attributes to ensure <code>Object[]</code> has consistent order (N elements).
 */
public class AttributeExtractor {

  /**
   */
  public static final class IncompatibleHeaderException extends Exception {
    /**
     * Construct an exception with given message.
     * @param message exception message.
     */
    public IncompatibleHeaderException(String message) {
      super(message);
    }
  }

  private final Annotation[] mAnnotations;
  private final Attribute[] mAttributes;

  /**
   * Constructs an attribute extractor for the given set of annotations and corresponding attributes.
   * @param annotations annotations responsible for extracting field values
   * @param attributes attributes responsible for translating field values into encoded representation
   */
  AttributeExtractor(Annotation[] annotations, Attribute[] attributes) {
    assert annotations.length == attributes.length;
    mAnnotations = annotations;
    mAttributes = attributes;
  }

  // Used during testing, assumes no need to re-use attributes across extractors
  AttributeExtractor(Annotation... annotations) {
    this(annotations, createAttributes(annotations));
  }

  static Annotation[] normalizeAnnotations(Collection<Annotation> annotations) {
    return normalizeAnnotations(annotations.toArray(new Annotation[0]));
  }

  static Annotation[] normalizeAnnotations(Annotation[] annotations) {
    assert annotations.length > 0;
    final SortedSet<Annotation> set = new TreeSet<>(Comparator.comparing(Annotation::getName));
    for (Annotation ann : annotations) {
      if (ann == null) {
        throw new NullPointerException("null annotation given");
      }
      if (set.contains(ann)) {
        throw new IllegalArgumentException("Duplicate annotation: " + ann.getName());
      }
      set.add(ann);
    }

    return set.toArray(new Annotation[0]);
  }

  static Attribute[] createAttributes(Annotation[] annotations) {
    final Attribute[] attributes = new Attribute[annotations.length];
    for (int i = 0; i < attributes.length; ++i) {
      attributes[i] = new Attribute(annotations[i].getName(), getMlDataType(annotations[i].getType()));
    }
    return attributes;
  }

  /**
   * Check that the annotations are compatible with the VCF header.
   *
   * @param header a VCF header
   * @throws IncompatibleHeaderException if the head and attributes don't match
   */
  public void checkHeader(VcfHeader header) throws IncompatibleHeaderException {
    final StringBuilder exceptionMessages = new StringBuilder();
    for (Annotation anno : mAnnotations) {
      final String message = anno.checkHeader(header);
      if (message != null) {
        exceptionMessages.append(message).append(LS);
      }
    }
    if (exceptionMessages.length() != 0) {
      throw new IncompatibleHeaderException(exceptionMessages.toString());
    }
  }

  /**
   * Returns empty dataset that complies with the attributes.
   * @return a dataset with attributes set up.
   */
  Dataset getDataset() {
    return new Dataset(mAttributes);
  }

  protected static MlDataType getMlDataType(AnnotationDataType dataType) {
    switch (dataType) {
      case BOOLEAN:
        return MlDataType.BOOLEAN;
      case INTEGER:
        return MlDataType.INTEGER;
      case DOUBLE:
        return MlDataType.DOUBLE;
      case STRING:
        return MlDataType.STRING;
      default:
        throw new IllegalArgumentException("Unrecognised Data Type: " + dataType);
    }
  }

  /**
   * Returns an instance object array for the given VCF record.
   * @param record a VCF record
   * @param sampleNumber the sample to extract value from.
   * @return array of attribute values
   */
  public double[] getInstance(VcfRecord record, int sampleNumber) {
    final double[] res = new double[mAnnotations.length];
    for (int i = 0; i < res.length; ++i) {
      try {
        res[i] = mAttributes[i].encodeValue(mAnnotations[i].getValue(record, sampleNumber));
      } catch (final NumberFormatException e) {
        throw new NoTalkbackSlimException("Problem parsing a number in a VCF record:\n" + record + "\n" + e);
      }
    }
    return res;
  }

  double[] getMissingValuesInstance() {
    final double[] nullInstance = new double[mAnnotations.length];
    Arrays.fill(nullInstance, Double.NaN);
    return nullInstance;
  }

  /**
   * Return a summary of the number of missing values.
   * @param data the input dataset
   * @return missing values summary
   */
  public String missingValuesReport(Dataset data) {
    final long[] counts = data.missingValueCounts();
    final int numAtts = data.getAttributes().length;
    final TextTable t = new TextTable(2, 2, TextTable.Align.RIGHT);
    t.setAlignment(TextTable.Align.LEFT, TextTable.Align.RIGHT);
    for (int k = 0; k < numAtts; ++k) {
      t.addRow(mAnnotations[k].getName(), String.valueOf(counts[k]));
    }
    return "Number of examples with missing values:" + LS + t;
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    for (Annotation att : mAnnotations) {
      if (sb.length() != 0) {
        sb.append(",");
      }
      sb.append(att.getName()).append("(").append(att.getType()).append(")");
    }
    return sb.toString();
  }

  /**
   * Creates an {@link AttributeExtractor} form the given input stream.
   * @param is an input stream
   * @return an {@link AttributeExtractor} from stream contents
   * @throws IOException if an error occurs reading data
   */
  public static AttributeExtractor load(InputStream is) throws IOException {
    final DataInputStream dis = new DataInputStream(is);
    final int numAnnotations = dis.readInt();
    final Annotation[] annotations = new Annotation[numAnnotations];
    for (int i = 0; i < numAnnotations; ++i) {
      annotations[i] = AnnotationLoader.load(dis);
    }
    return new AttributeExtractor(annotations, createAttributes(annotations));
  }

  /**
   * Saves this {@link AttributeExtractor} to the given output stream.
   * @param os output stream to write to
   * @throws IOException if an error occurs writing data
   */
  public void save(OutputStream os) throws IOException {
    final DataOutputStream dos = new DataOutputStream(os);
    dos.writeInt(mAnnotations.length);
    try {
      for (Annotation anno : mAnnotations) {
        anno.save(dos);
      }
    } finally {
      dos.flush();
    }
  }

  protected static AnnotationDataType getCompatibleType(MetaType mt) {
    switch (mt) {
      case INTEGER:
        return AnnotationDataType.INTEGER;
      case FLOAT:
        return AnnotationDataType.DOUBLE;
      case STRING:
      case CHARACTER:
        return AnnotationDataType.STRING;
      case FLAG:
        return AnnotationDataType.BOOLEAN;
      default:
        return null;
    }
  }

}
