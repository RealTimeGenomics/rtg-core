/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.variant.avr;

import java.io.File;
import java.io.IOException;
import java.io.OutputStreamWriter;

import com.rtg.ml.Attribute;
import com.rtg.ml.Dataset;
import com.rtg.ml.Instance;
import com.rtg.util.io.BaseFile;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LineWriter;

/**
 * Writes training data in ARFF format
 */
public class ArffDatasetWriter extends MlAvrModelBuilder {

  // String attributes with fewer than this many distinct values will be written as nominal attributes
  static final int MAX_NOMINAL = 20;

  private Dataset mDataset = null;

  /**
   * Create a Machine Learning AVR model by extracting values from the given attribute fields.
   * @param formatAttributes FORMAT field attributes to use
   * @param infoAttributes INFO field attributes to use
   * @param derivedAttributes derived field attributes to use
   */
  public ArffDatasetWriter(String[] formatAttributes, String[] infoAttributes, String[] derivedAttributes) {
    super(formatAttributes, infoAttributes, derivedAttributes);
  }


  @Override
  public void build(VcfDataset... vcfDatasets) throws IOException {
    scanVcfHeaders(vcfDatasets);
    final Annotation[] annotations = createAnnotations();
    final Attribute[] attributes = AttributeExtractor.createAttributes(annotations);
    mDataset = extractDataset(attributes, vcfDatasets);
  }

  /**
   * Instead of outputting a model, outputs the training data in ARFF format
   */
  @Override
  public void save(File arffFile) throws IOException {
    final Annotation[] annotations = createAnnotations();
    final BaseFile bf = FileUtils.getBaseFile(arffFile, false, ".arff");
    final Attribute[] attributes = mDataset.getAttributes();
    try (final LineWriter w = new LineWriter(new OutputStreamWriter(FileUtils.createOutputStream(bf, "")))) {
      w.writeln("@relation '" + bf.getBaseFile().getName() + "'");
      w.writeln("");

      for (int i = 0; i < annotations.length; i++) {
        w.writeln(getArffAttributeText(annotations[i], attributes[i]));
      }
      w.writeln("@attribute class {true,false}");

      w.writeln("");
      w.writeln("@data");
      w.writeln("");

      for (Instance theInstance : mDataset.getInstances()) {
        final double[] instance = theInstance.instance();
        for (int i = 0; i < instance.length; ++i) {
          if (i != 0) {
            w.write(",");
          }
          w.write(getValueAsArffString(attributes[i].decodeValue(instance[i])));
        }
        w.write(",");
        w.write(Boolean.valueOf(theInstance.isPositive()).toString());
        if (theInstance.weight() != 1.0) {
          w.write(",{");
          w.write(Double.toString(theInstance.weight()));
          w.write("}");
        }
        w.writeln("");
      }
      w.flush();
    }
  }

  protected static String getArffAttributeText(Annotation annotation, Attribute att) {
    final StringBuilder sb = new StringBuilder("@attribute ");
    sb.append(annotation.getName());
    switch (annotation.getType()) {
      case INTEGER:
      case DOUBLE:
        sb.append(" numeric");
        break;
      case BOOLEAN:
        sb.append(" {").append(Boolean.TRUE).append(",").append(Boolean.FALSE).append("}");
        break;
      case STRING:
        if (att.nominalSize() < MAX_NOMINAL) {
          sb.append(" {");
          for (int i = 0; i < att.nominalSize(); i++) {
            if (i != 0) {
              sb.append(",");
            }
            sb.append(att.decodeValue(i));
          }
          sb.append("}");
        } else {
          sb.append(" string");
        }
        break;
      default:
        throw new UnsupportedOperationException();
    }
    return sb.toString();
  }

  protected static String getValueAsArffString(Object value) {
    return value == null ? "?" : value.toString(); //.replace(',', '_');
  }
}
