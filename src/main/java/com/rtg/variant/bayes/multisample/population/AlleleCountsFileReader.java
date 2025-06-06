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

package com.rtg.variant.bayes.multisample.population;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.TabixLineReader;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.vcf.VariantType;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;
import com.rtg.vcf.header.InfoField;
import com.rtg.vcf.header.VcfHeader;

/**
 * Reads either a <code>VCF</code> or an allele counts file into a string array ready for parsing
 * Warnings handling not thread safe.
 */
@TestClass(value = {"com.rtg.variant.bayes.multisample.population.AlleleCountsFileReaderTest", "com.rtg.variant.bayes.multisample.population.PopulationHwHypothesesCreatorTest"})
public final class AlleleCountsFileReader implements Closeable {

  private static final String AC_INFO_ID = "AC";
  private static final String AN_INFO_ID = "AN";

  private final TabixLineReader mReader;
  private final VcfReader mVcfReader;
  private final boolean mSupportsANAC;

  private AlleleCounts mCurrent;
  private String mCurrentReference;

  private Integer mWarnings = 0;

  /**
   * @param reader tabix line reader for non-VCF input
   */
  private AlleleCountsFileReader(TabixLineReader reader) {
    mReader = reader;
    mVcfReader = null;
    mSupportsANAC = false;
  }

  /**
   * @param vcfReader reader for VCF input
   */
  private AlleleCountsFileReader(VcfReader vcfReader) {
    mReader = null;
    mVcfReader = vcfReader;
    mSupportsANAC = checkSupportsACandAN(mVcfReader.getHeader());
  }

  /**
   * Create an AlleleCountsFileReader
   * @param f the allele counts file to read from
   * @param ranges the ranges restriction to apply, if any
   * @return an AlleleCountsFileReader to read from
   * @throws IOException if bad things occur
   */
  public static AlleleCountsFileReader openAlleleCountReader(File f, ReferenceRanges<String> ranges) throws IOException {
    if (VcfUtils.isVcfExtension(f)) {
      return new AlleleCountsFileReader(VcfReader.openVcfReader(f, ranges));
    } else {
      return new AlleleCountsFileReader(new TabixLineReader(f, TabixIndexer.indexFileName(f), ranges));
    }
  }


  public AlleleCounts getCurrent() {
    return mCurrent;
  }

  public String getCurrentReference() {
    return mCurrentReference;
  }

  /**
   * Consume the next record if possible
   * @return true if successful, or false in no further records available
   * @throws IOException if bad things occur
   */
  public boolean next() throws IOException {
    if (mVcfReader != null) { //entirely certain this vcf section can made a lot nicer but don't have time right now
      while (mVcfReader.hasNext()) {
        final VcfRecord rec = mVcfReader.next();
        mCurrent = vcfRecordToAlleleCountLine(rec);
        if (mCurrent != null) {
          mCurrentReference = rec.getSequenceName();
          return true;
        }
      }
    } else { //this is the 'cut down text file' version
      String line;
      while ((line = mReader.readLine()) != null) {
        mCurrent = parseTextFileFormat(line);
        if (mCurrent != null) {
          //current reference is set during parse
          return true;
        }
      }
    }
    mCurrent = null;
    mCurrentReference = null;
    return false;
  }

  private AlleleCounts parseTextFileFormat(String line) {
    final String[] splitLine = StringUtils.split(line, '\t');
    if (splitLine.length < 4) {
      throw new NoTalkbackSlimException("Line must contain at least <templateId> <position> <ref allele> <allelecount>: " + line);
    }
    final String refCall = splitLine[2];
    if (VcfRecord.MISSING.equals(refCall)) {
      Diagnostic.warning("Empty value in reference field not supported: " + line);
      return null;
    }

    final Map<String, Integer> countsMap = new HashMap<>(splitLine.length);
    for (int i = 2; i < splitLine.length; i += 2) {
      final String allele = splitLine[i];
      final Integer count;
      try {
        count = Integer.valueOf(splitLine[i + 1]);
      } catch (final NumberFormatException nfe) {
        throw new NoTalkbackSlimException(splitLine[i + 1] + " not an integer in line: " + line);
      }
      countsMap.put(allele, count);
    }

    final Integer position;
    try {
      position = Integer.valueOf(splitLine[1]);
    } catch (final NumberFormatException nfe) {
      throw new NoTalkbackSlimException("Position: " + splitLine[1] + " not an integer in line: " + line);
    }
    mCurrentReference = splitLine[0];
    return new AlleleCounts(position - 1, countsMap, refCall);
  }

  /**
   * Checks if the <code>VCF</code> header in a <code>VCF</code> reader contains info field declarations for both AC and AN
   * @param vcfHeader the <code>VCF</code> header to check
   * @return true if both AN and AC are defined info fields in the header
   */
  public static boolean checkSupportsACandAN(VcfHeader vcfHeader) {
    final List<InfoField> infoLines = vcfHeader.getInfoLines();
    boolean supportsAN = false;
    boolean supportsAC = false;
    for (final InfoField infoField : infoLines) {
      if (AN_INFO_ID.equals(infoField.getId())) {
        supportsAN = true;
      } else if (AC_INFO_ID.equals(infoField.getId())) {
        supportsAC = true;
      }
      if (supportsAC && supportsAN) {
        break;
      }
    }
    return supportsAC && supportsAN;
  }

  /**
   * Converts a <code>VCF</code> record into a string conforming to the allele count file format
   * @param vcfRecord a <code>VCF</code> record
   * @return a string representation of the <code>VCF</code> record in the allele count file format, or null if not possible
   */
  public AlleleCounts vcfRecordToAlleleCountLine(VcfRecord vcfRecord) {
    if (VcfRecord.MISSING.equals(vcfRecord.getRefCall())) {
      ++mWarnings;
      if (mWarnings < 10) {
        Diagnostic.warning("Empty value in reference field not supported: " + vcfRecord);
      }
      return null;
    } else if (vcfRecord.isFiltered()) {
      return null;
    }

    for (final String alt : vcfRecord.getAltCalls()) {
      if (VariantType.getType(vcfRecord.getRefCall(), alt).isSvType()) {
        return null;
      }
    }
    final Map<String, Integer> countsMap;
    if (mSupportsANAC) {
      if (vcfRecord.hasInfo(AN_INFO_ID) && vcfRecord.hasInfo(AC_INFO_ID)) {
        final String[] anValue = vcfRecord.getInfoSplit(AN_INFO_ID);
        final String[] acValue = vcfRecord.getInfoSplit(AC_INFO_ID);

        if (anValue.length != 1) {
          ++mWarnings;
          if (mWarnings < 10) {
            Diagnostic.warning("INFO field " + AN_INFO_ID + " contains too many values: " + vcfRecord);
          }
          return null;
        }
        final Integer totalAlleleNumber;
        try {
          totalAlleleNumber = Integer.valueOf(anValue[0]);
        } catch (final NumberFormatException nfe) {
          ++mWarnings;
          if (mWarnings < 10) {
            Diagnostic.warning("INFO field " + AN_INFO_ID + " not an integer: " + vcfRecord);
          }
          return null;
        }
        countsMap = new HashMap<>();
        int totalAllelesSeen = 0;

        if (acValue.length != vcfRecord.getAltCalls().size()) {
          Diagnostic.warning("Allele count field length not equal to number of alternate alleles, was " + acValue.length + ", expected " + vcfRecord.getAltCalls().size());
          return null;
        }
        for (int i = 0; i < acValue.length; ++i) {
          try {
            final Integer thisAlleleCount = Integer.valueOf(acValue[i]);
            countsMap.put(vcfRecord.getAltCalls().get(i), thisAlleleCount);
            totalAllelesSeen += thisAlleleCount;
          } catch (final NumberFormatException nfe) {
            ++mWarnings;
            if (mWarnings < 10) {
              Diagnostic.warning("INFO field " + AC_INFO_ID + " contained a non-integer: " + vcfRecord);
            }
            return null;
          }
        }
        countsMap.put(vcfRecord.getRefCall(), totalAlleleNumber - totalAllelesSeen);

      } else {
        countsMap = initCountsTheHardWay(vcfRecord);
      }
    } else {
      countsMap = initCountsTheHardWay(vcfRecord);
    }
    return new AlleleCounts(vcfRecord.getStart(), countsMap, vcfRecord.getRefCall());
  }

  private Map<String, Integer> initCountsTheHardWay(VcfRecord vcfRecord) {
    final Map<String, Integer> countsMap = new HashMap<>();

    final int[] counts = new int[vcfRecord.getAltCalls().size() + 1]; // +1 for the reference allele.
    for (int i = 0; i < vcfRecord.getNumberOfSamples(); ++i) {
      final String gt = vcfRecord.getFormat(VcfUtils.FORMAT_GENOTYPE).get(i);
      for (final int gti : VcfUtils.splitGt(gt)) {
        if (gti != -1) { // to ignore . alleles in GT like ./1 which is valid (e.g. on sex chromosome PAR regions, depending on representation).
          if (gti >= counts.length) {
            throw new NoTalkbackSlimException("GT field referenced ALT which does not exist on line: " + vcfRecord);
          }
          counts[gti]++;
        }
      }
    }

    countsMap.put(vcfRecord.getRefCall(), counts[0]);
    for (int i = 0; i < vcfRecord.getAltCalls().size(); ++i) {
      final String allele = vcfRecord.getAltCalls().get(i);
      if (countsMap.containsKey(allele)) {
        duplicateAlleleWarning(allele, vcfRecord);
        continue;
      }
      countsMap.put(vcfRecord.getAltCalls().get(i), counts[i + 1]); //i-th alt = i-th + 1 count (0 is for ref)
    }

    return countsMap;
  }

  private synchronized void duplicateAlleleWarning(String allele, VcfRecord vcfRecord) {
    ++mWarnings;
    if (mWarnings < 10) {
      Diagnostic.warning("Ignoring duplicate allele: " + allele + " in VCF record: " + vcfRecord);
    }
  }

  @Override
  public void close() throws IOException {
    if (mVcfReader != null) {
      mVcfReader.close();
    }
    if (mReader != null) {
      mReader.close();
    }
    if (mWarnings > 0) {
      Diagnostic.warning(mWarnings + " total warnings during allele counts reading.");
    }
  }
}
