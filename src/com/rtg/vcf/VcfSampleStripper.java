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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.vcf.header.VcfHeader;

/**
 * Removes unwanted samples from a VCF record
 */
public class VcfSampleStripper implements VcfAnnotator {

  private final boolean mRemoveAll;
  private final boolean mKeepMode;
  private final Set<String> mSamples;
  private int[] mSampleIdsToRemove = null;

  /**
   * Remove all samples from header and records
   * @param removeAll false if you don't actually want to do it for some reason.
   */
  VcfSampleStripper(boolean removeAll) {
    mRemoveAll = removeAll;
    mKeepMode = false;
    mSamples = null;
  }

  /**
   * Keep or remove a selected set of samples from header and records
   * @param sampleList the list of sample names
   * @param keep true to keep values in the list, false to remove them
   */
  VcfSampleStripper(Set<String> sampleList, boolean keep) {
    mRemoveAll = false;
    mKeepMode = keep;
    mSamples = sampleList;
  }

  @Override
  public void updateHeader(VcfHeader header) {
    if (mRemoveAll) {
      header.removeAllSamples();
    } else {
      final HashSet<String> samplesToRemove = new HashSet<>();
      for (String sample : header.getSampleNames()) {
        if (mKeepMode ^ mSamples.contains(sample)) {
          samplesToRemove.add(sample);
        }
      }

      mSampleIdsToRemove = new int[samplesToRemove.size()];
      int i = 0;
      for (String sample : samplesToRemove) {
        final Integer sampleId = header.getSampleIndex(sample);
        if (sampleId == null) {
          throw new NoTalkbackSlimException("Could not find sample name: " + sample + " in VCF header");
        }
        mSampleIdsToRemove[i] = sampleId;
        i++;
      }
      Arrays.sort(mSampleIdsToRemove); //this so that when we come to remove from the format sample lists,

      header.removeSamples(samplesToRemove);
    }
    if (header.getNumberOfSamples() == 0) {
      header.getFormatLines().clear();
    }
  }

  @Override
  public void annotate(VcfRecord rec) {
    if (mRemoveAll) {
      rec.getFormatAndSample().clear();
      rec.setNumberOfSamples(0);
      return;
    } else if (mSamples == null || mSamples.size() == 0) {
      return;
    }

    if (mSampleIdsToRemove == null) {
      throw new RuntimeException("Call updateHeader first.");
    }

    int newNumSamples = rec.getNumberOfSamples();
    boolean first = true;
    for (Map.Entry<String, ArrayList<String>> formatValue : rec.getFormatAndSample().entrySet()) {
      //remove each sample from each format
      for (int j = mSampleIdsToRemove.length - 1; j >= 0; j--) { //backwards to avoid changing index values as we remove items
        formatValue.getValue().remove(mSampleIdsToRemove[j]);
      }
      if (first) {
        newNumSamples = formatValue.getValue().size();
        first = false;
      }
    }
    rec.setNumberOfSamples(newNumSamples);
    if (newNumSamples == 0) {
      rec.getFormatAndSample().clear();
    }
  }
}
