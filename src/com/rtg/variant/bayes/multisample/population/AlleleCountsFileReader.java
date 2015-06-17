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

package com.rtg.variant.bayes.multisample.population;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
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
    final Map<String, Integer> countsMap = new HashMap<>();

    for (int i = 2; i < splitLine.length; i += 2) {
      final String allele = splitLine[i];
      final Integer count;
      try {
        count = Integer.parseInt(splitLine[i + 1]);
      } catch (final NumberFormatException nfe) {
        throw new NoTalkbackSlimException(splitLine[i + 1] + " not an integer in line: " + line);
      }
      countsMap.put(allele, count);
    }

    final Integer position;
    try {
      position = Integer.parseInt(splitLine[1]);
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
      if (infoField.getId().equals(AN_INFO_ID)) {
        supportsAN = true;
      } else if (infoField.getId().equals(AC_INFO_ID)) {
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
      mWarnings++;
      if (mWarnings < 10) {
        Diagnostic.warning("Empty value in reference field not supported: " + vcfRecord.toString());
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
      final Map<String, ArrayList<String>> info = vcfRecord.getInfo();
      if (info.containsKey(AN_INFO_ID) && info.containsKey(AC_INFO_ID)) {
        final ArrayList<String> anValue = info.get(AN_INFO_ID);
        final ArrayList<String> acValue = info.get(AC_INFO_ID);

        if (anValue.size() != 1) {
          mWarnings++;
          if (mWarnings < 10) {
            Diagnostic.warning("INFO field " + AN_INFO_ID + " contains too many values: " + vcfRecord.toString());
          }
          return null;
        }
        final Integer totalAlleleNumber;
        try {
          totalAlleleNumber = Integer.parseInt(anValue.get(0));
        } catch (final NumberFormatException nfe) {
          mWarnings++;
          if (mWarnings < 10) {
            Diagnostic.warning("INFO field " + AN_INFO_ID + " not an integer: " + vcfRecord.toString());
          }
          return null;
        }
        countsMap = new HashMap<>();
        int totalAllelesSeen = 0;

        if (acValue.size() != vcfRecord.getAltCalls().size()) {
          Diagnostic.warning("Allele count field length not equal to number of alternate alleles, was " + acValue.size() + ", expected " + vcfRecord.getAltCalls().size());
          return null;
        }
        for (int i = 0; i < acValue.size(); i++) {
          try {
            final Integer thisAlleleCount = Integer.parseInt(acValue.get(i));
            countsMap.put(vcfRecord.getAltCalls().get(i), thisAlleleCount);
            totalAllelesSeen += thisAlleleCount;
          } catch (final NumberFormatException nfe) {
            mWarnings++;
            if (mWarnings < 10) {
              Diagnostic.warning("INFO field " + AC_INFO_ID + " contained a non-integer: " + vcfRecord.toString());
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
    for (int i = 0; i < vcfRecord.getNumberOfSamples(); i++) {
      final String gt = vcfRecord.getFormatAndSample().get(VcfUtils.FORMAT_GENOTYPE).get(i);
      for (final int gti : VcfUtils.splitGt(gt)) {
        if (gti != -1) { // to ignore . alleles in GT like ./1 which is valid (e.g. on sex chromosome PAR regions, depending on representation).
          if (gti >= counts.length) {
            throw new NoTalkbackSlimException("GT field referenced ALT which does not exist on line: " + vcfRecord.toString());
          }
          counts[gti]++;
        }
      }
    }

    countsMap.put(vcfRecord.getRefCall(), counts[0]);
    for (int i = 0; i < vcfRecord.getAltCalls().size(); i++) {
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
    mWarnings++;
    if (mWarnings < 10) {
      Diagnostic.warning("Ignoring duplicate allele: " + allele + " in VCF record: " + vcfRecord.toString());
    }
  }

  /**
   * Close the underlying reader
   * @throws IOException if bad things occur
   */
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
