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
