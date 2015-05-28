package com.rtg.vcf;

/**
 * @author Sean A. Irvine
 */
public interface VcfFilter {
  boolean accept(VcfRecord record);
}
