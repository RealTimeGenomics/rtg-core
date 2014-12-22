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
package com.rtg.vcf.header;

/**
 * Interface for <code>VCF</code> header (meta) lines using id fields
 * @param <T> sub class type
 */
public interface IdField<T extends IdField<T>> {

  /**
   * @return the id of this field
   */
  String getId();

  /**
   * If this and the other field can be represented by a common <code>IdField</code> then
   * return that field. Otherwise return null to indicate they are not compatible.
   * @param other the field to compare to
   * @return the field that equivalently represents both this and the other field. null if fields are no compatible
   */
  T superSet(T other);

  /**
   * @return in file string representation of the meta line
   */
  String toString();

}
