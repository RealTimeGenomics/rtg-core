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
package com.rtg.vcf;

import com.rtg.vcf.header.VcfHeader;

/**
 */
public interface VcfFilter {

  /**
   * Process a VCF record and record statistics.
   * @param record the VCF record to filter
   * @return false if the record is unacceptable and should be filtered
   */
  boolean accept(VcfRecord record);

  /**
   * Set the VCF header to be used with this filter.
   * @param header VCF header
   */
  void setHeader(VcfHeader header);


}
