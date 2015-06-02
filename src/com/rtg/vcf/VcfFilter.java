package com.rtg.vcf;

import com.rtg.vcf.header.VcfHeader;

/**
 * @author Sean A. Irvine
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
