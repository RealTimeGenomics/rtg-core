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

package com.rtg.variant.avr;

import java.io.DataOutputStream;
import java.io.IOException;

import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.header.VcfHeader;

/**
 * VCF annotation specification and extraction interface.
 */
public interface Annotation {

   /**
   * Gets the name of the annotation.
   * @return name of annotation
   */
  String getName();

  /**
   * Returns the type of the annotation.
   * @return annotation type.
   */
  AnnotationDataType getType();

  /**
   * Returns the value for an annotation from a given {@link VcfRecord}.
   * For sample specific annotations, the value in the first sample is returned.
   * @param record record to extract value from.
   * @param sampleNumber the number of the sample to extract the value from.
   * @return value for this annotation from the given record, or null if no value could be calculated.
   */
  Object getValue(VcfRecord record, int sampleNumber);

  /**
   * Ensures any required annotations are declared in the given header.
   * @param header a VCF header
   * @return null if all good, other wise an error message.
   */
  String checkHeader(VcfHeader header);

  /**
   * Save the annotation to given data output stream.
   * @param dos a data output stream
   * @throws IOException if an error occurs saving
   */
  void save(DataOutputStream dos) throws IOException;
}
