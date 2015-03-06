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
package com.rtg.sam;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.SortedMap;
import java.util.TreeMap;

import com.rtg.mode.DNA;
import com.rtg.mode.DnaUtils;
import com.rtg.ngs.NgsParams;
import com.rtg.reader.PrereadType;
import com.rtg.reader.ReadHelper;
import com.rtg.reader.SdfId;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.io.FileUtils;

import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFormatException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

/**
 * Verifies a SAM file produced from paired end is sane.
 */
public final class SamValidator {

  private final boolean mValidate;
  private final boolean mShowConsensus;
  private final boolean mPrintHistograms;
  private final PrintStream mOut;
  private final PrintStream mErr;
  private SequencesReader mLeftReader;
  private SequencesReader mRightReader;
  private final boolean mIgnoreFragmentSizeProblems;
  private final boolean mPerFileStats;

  private final SuperCigarValidator mSuperCigarValidator;
  private final HashSet<String> mExpectedMates = new HashSet<>();
  private long mFormatProblems = 0;
  private long mExceptionProblems = 0;
  private boolean mSingleEnded = false;

  private final SamStatsVariables mTotalVariables = new SamStatsVariables();
  private SamStatsVariables mCurrentVariables = new SamStatsVariables();

  private final int mMismatchPenalty;
  private final int mUnknownsPenalty;
  private final int mGapOpenPenalty;
  private final int mGapExtendPenalty;
  private final boolean mPenaltiesSet;

  //CHECKSTYLE:OFF
  SamValidator(PrintStream out, PrintStream err, boolean validate, boolean showConsensus, boolean printHistograms, boolean ignoreCgFragmentSize, boolean perFileStats, NgsParams params, boolean penaltiesSet) {
  //CHECKSTYLE:ON
    mOut = out;
    mErr = err;
    mValidate = validate;
    mShowConsensus = showConsensus;
    mPrintHistograms = printHistograms;
    mIgnoreFragmentSizeProblems = ignoreCgFragmentSize;
    mPerFileStats = perFileStats;

    mPenaltiesSet = penaltiesSet;
    mMismatchPenalty = params.substitutionPenalty();
    mUnknownsPenalty = params.unknownsPenalty();
    mGapOpenPenalty = params.gapOpenPenalty();
    mGapExtendPenalty = params.gapExtendPenalty();

    mSuperCigarValidator = new SuperCigarValidator(mUnknownsPenalty > 0 ? 1 : 0);
  }

  /**
   * Construct a Sam Validator
   * @param out the output stream
   * @param err the error stream
   * @param validate true to validate each record
   * @param showConsensus true to display consensus information
   * @param printHistograms true to print histograms
   * @param perFileStats true to output per file statistics
   * @param params {@link com.rtg.ngs.NgsParams} for current run
   * @param penaltiesSet true if the user has specified the penalties (for alignment score checking)
   */
  public SamValidator(PrintStream out, PrintStream err, boolean validate, boolean showConsensus, boolean printHistograms, boolean perFileStats, NgsParams params, boolean penaltiesSet) {
    this(out, err, validate, showConsensus, printHistograms, false, perFileStats, params, penaltiesSet);
  }

  void checkSAMAlign(File templateDir, Collection<File> samFiles, File leftReadsDir, File rightReadsDir) throws IOException {
    try {
      try (SequencesReader templateReader = SequencesReaderFactory.createDefaultSequencesReaderCheckEmpty(templateDir)) {
        boolean cgData = false;
        SdfId readsGuid = new SdfId(0);
        try {
          mLeftReader = SequencesReaderFactory.createMemorySequencesReader(leftReadsDir, false, LongRange.NONE);
        } catch (final FileNotFoundException e) {
          throw new NoTalkbackSlimException(ErrorType.SDF_INDEX_NOT_VALID, leftReadsDir.toString());
        }
        if (mLeftReader != null) {
          cgData = mLeftReader.getPrereadType() == PrereadType.CG;
          readsGuid = mLeftReader.getSdfId();
        }
        try {
          mRightReader = SequencesReaderFactory.createMemorySequencesReader(rightReadsDir, false, LongRange.NONE);
        } catch (final FileNotFoundException e) {
          throw new NoTalkbackSlimException(ErrorType.SDF_INDEX_NOT_VALID, rightReadsDir.toString());
        }
        if (mRightReader != null) {
          if (!mRightReader.getSdfId().check(readsGuid)) {
            throw new NoTalkbackSlimException(ErrorType.FILE_READ_ERROR, "Left and right reads have different GUIDs - are from different runs.");
          }
          if (cgData && mRightReader.getPrereadType() != PrereadType.CG) {
            throw new NoTalkbackSlimException(ErrorType.FILE_READ_ERROR, "Left reads are CG data, but right reads are not.");
          }
        }
        if (mLeftReader != null && mRightReader == null) {
          mSingleEnded = true;
        }
        final int[] countPerRead = mLeftReader == null ? mRightReader == null ? null : new int[(int) mRightReader.numberSequences()] : new int[(int) mLeftReader.numberSequences()];
        final boolean[] pairedRead = mLeftReader == null ? mRightReader == null ? null : new boolean[(int) mRightReader.numberSequences()] : new boolean[(int) mLeftReader.numberSequences()];
        for (final File samFile : samFiles) {
          templateReader.seek(0);
          mExpectedMates.clear();
          try (InputStream bis = FileUtils.createInputStream(samFile, false)) {
            try (SamReader read = new SAMFileReader(bis)) {
              processRecords(templateReader, read, cgData, countPerRead, pairedRead, mCurrentVariables);
              if (mValidate) {
                for (final String mate : mExpectedMates) {
                  mErr.println("Missing mate: " + mate);
                }
              }

              if (mPerFileStats || samFiles.size() == 1) {
                printStats(samFile, mCurrentVariables);
              }
            }
          }
          mTotalVariables.addToTotal(mCurrentVariables);
          mCurrentVariables = new SamStatsVariables();
        }
        if (samFiles.size() > 1) {
          printStats(null, mTotalVariables);
        }
        printStats2(countPerRead, pairedRead, mOut, mSingleEnded);
      } finally {
        if (mLeftReader != null) {
          mLeftReader.close();
        }
        if (mRightReader != null) {
          mRightReader.close();
        }
      }
    } catch (final FileNotFoundException e) {
      throw new NoTalkbackSlimException(ErrorType.SDF_INDEX_NOT_VALID, templateDir.toString());
    }
  }

  private void printStats(File samFile, SamStatsVariables variables) {
    if (samFile != null) {
      mOut.println("Stats for file: " + samFile.getPath());
    } else {
      mOut.println("Overall Statistics");
    }
    printBasicStats(variables.mTotalRecords, variables.mUnmappedRecords, mFormatProblems);
    printStats(variables);
    if (mPrintHistograms) {
      printDistributions(variables, mOut);
    }
    mOut.println();
  }

  private void validate(SAMRecord samRec, int prevTemplatePosition, int readId, boolean first, boolean cgData) {
    if (samRec.getAlignmentStart() < prevTemplatePosition) {
      mErr.println("Ordering error: " + samRec.getSAMString().trim());
    }
    if (samRec.getReadPairedFlag()) {
      if (samRec.getFirstOfPairFlag() == samRec.getSecondOfPairFlag()) {
        mErr.println("Mates from same side " + samRec.getSAMString().trim());
      }
      final String mateRef = samRec.getMateReferenceName();
      if (!samRec.getMateUnmappedFlag() && !"*".equals(mateRef)) {
        if (!mateRef.equals(samRec.getReferenceName())) {
          mErr.println("Mate ref name not same as record ref name " + samRec.getSAMString().trim());
        }
        final String mateKey = samRec.getMateAlignmentStart() + ":" + samRec.getAlignmentStart() + ":" + samRec.getReadName() + samRec.getSecondOfPairFlag() + samRec.getReadNegativeStrandFlag() + samRec.getMateNegativeStrandFlag() + samRec.getInferredInsertSize();
        final String key = samRec.getAlignmentStart() + ":" + samRec.getMateAlignmentStart() + ":" + samRec.getReadName() + samRec.getFirstOfPairFlag() + samRec.getMateNegativeStrandFlag() + samRec.getReadNegativeStrandFlag() + -samRec.getInferredInsertSize();
        if (!mExpectedMates.remove(mateKey)) {
          mExpectedMates.add(key);
        }
      }   //else the mate is unmapped...
    }
    if (!matchesRawRead(read(readId, first), qual(readId, first), samRec, cgData, mIgnoreFragmentSizeProblems)) {
      mErr.println("Read doesn't match expected value from SDF " + samRec.getSAMString().trim());
    }
    final Integer ih = samRec.getIntegerAttribute(SamUtils.ATTRIBUTE_IH);
    if (ih != null && ih <= 0) {
      mErr.println("IH value invalid " + samRec.getSAMString().trim());
    }
    final Integer nh = samRec.getIntegerAttribute(SamUtils.ATTRIBUTE_NH);
    if (nh != null) {
      if (nh <= 0) {
        mErr.println("NH value invalid " + samRec.getSAMString().trim());
      } else if (ih != null && nh < ih) {
        mErr.println("NH should be greater than or equal to IH " + samRec.getSAMString().trim());
      }
    }
  }

  int[] processRecords(SequencesReader templateReader, SamReader read, boolean cgData, int[] countPerRead, boolean[] pairedRead, SamStatsVariables variables) throws IOException {
    int prevTemplatePosition = 0;
    String prevTemplateName = null;
    byte[] currTemplate = null;
    PileUp pileUp = null;
    boolean firstRecord = true;

    for (final Iterator<SAMRecord> it = read.iterator(); it.hasNext(); ) {
      final SAMRecord samRec;
      try {
        samRec = it.next();
      } catch (final SAMFormatException e) {
        if (++mFormatProblems < 10) {
          mErr.println(e.getMessage());
        }
        continue;
      } catch (final RuntimeException e) {
        mExceptionProblems++;
        if (mExceptionProblems < 100) {
          mErr.println(e.getMessage());
          continue;
        } else {
          throw e;
        }
      }
      if (firstRecord) {
        if (samRec.getHeader().getSortOrder() == SortOrder.unsorted) {
          throw new NoTalkbackSlimException(ErrorType.FILE_READ_ERROR, "SAM file must be sorted.");
        }
        firstRecord = false;
      }
      mCurrentVariables.mTotalRecords++;
      final String refName = samRec.getReferenceName();
      final int readId;
      try {
        readId = (mRightReader != null || mLeftReader != null) ? Integer.parseInt(samRec.getReadName()) : -1;
      } catch (final NumberFormatException nfe) {
        throw new NoTalkbackSlimException(ErrorType.FILE_READ_ERROR, "Only unrenamed sam files are supported.");
      }
      if (!"*".equals(samRec.getBaseQualityString())) {
        if (samRec.getBaseQualities().length != samRec.getReadLength()) {
          mErr.print("Read length and quality length differ " + samRec.getSAMString().trim());
        }
      }
      final boolean first = !samRec.getReadPairedFlag() || samRec.getFirstOfPairFlag();
      if (mValidate) {
        if (first && mLeftReader != null && readId >= mLeftReader.numberSequences() || !first && mRightReader != null && readId >= mRightReader.numberSequences()) {
          throw new NoTalkbackSlimException("Reads SDF doesn't match sam file - not enough reads");
        }
      }
      if ("*".equals(refName)) {        // unmapped record
        mCurrentVariables.mUnmappedRecords++;
        if (mValidate) {
          if (!matchesRawRead(read(readId, first), qual(readId, first), samRec, cgData, mIgnoreFragmentSizeProblems)) {
            mErr.println("Read doesn't match expected value from SDF file " + samRec.getSAMString().trim());
          }
        }
        // TODO probably more stuff can be checked here, and may need to handle unmated stuff
      } else {  // mapped record
        if (!refName.equals(prevTemplateName)) {
          prevTemplateName = refName;
          prevTemplatePosition = 0;
          if (currTemplate != null) {
            accumulatePileUp(pileUp, currTemplate.length);
          }
          pileUp = null;
          while (!templateReader.currentName().equals(samRec.getReferenceName())) {
            if (!templateReader.nextSequence()) {
              throw new NoTalkbackSlimException(ErrorType.FILE_READ_ERROR, "Sequence for " + samRec.getReferenceName() + " not found in template.");
            }
            if (!templateReader.currentName().equals(samRec.getReferenceName())) {
              mErr.println("Template \"" + templateReader.currentName() + "\" had no hits");
            }
          }
          currTemplate = new byte[templateReader.currentLength()];
          mSuperCigarValidator.setTemplate(currTemplate, currTemplate.length);
          if (mShowConsensus) {
            pileUp = new PileUp(currTemplate.length);
          }
          templateReader.readCurrent(currTemplate);
        }
        final int expectedRet = isAtExpectedRef(currTemplate, samRec, pileUp);
        if (expectedRet == -1) {
          mErr.println("Alignment mismatch " + samRec.getSAMString().trim());
        }
        if (samRec.getReadPairedFlag() && (samRec.getAlignmentStart() < samRec.getMateAlignmentStart() || (samRec.getAlignmentStart() == samRec.getMateAlignmentStart() && samRec.getFirstOfPairFlag()))) {
          mCurrentVariables.mPairOrientations[(samRec.getFlags() >> 4) & 3]++;
        }
        if (mValidate) {
          validate(samRec, prevTemplatePosition, readId, first, cgData);
          if (mPenaltiesSet) {
            checkAlignmentScore(samRec, variables.mAlignmentScores, expectedRet);
          }
        }
        accumulateAlignmentCounts(samRec, variables.mIH, SamUtils.ATTRIBUTE_IH);
        accumulateAlignmentCounts(samRec, variables.mNH, SamUtils.ATTRIBUTE_NH);
        if (countPerRead != null) {
          countPerRead[readId]++;
        }
        if (pairedRead != null && samRec.getReadPairedFlag() && samRec.getProperPairFlag() && !samRec.getMateUnmappedFlag()) {
          pairedRead[readId] = true;
        }
        final int readLength = samRec.getReadLength();
        if (mCurrentVariables.mMaxReadLength < readLength) {
          mCurrentVariables.mMaxReadLength = readLength;
        }
        if (mCurrentVariables.mMinReadLength > readLength) {
          mCurrentVariables.mMinReadLength = readLength;
        }
        if (samRec.getReadPairedFlag() && samRec.getProperPairFlag()) {
          final Integer insertSize = samRec.getInferredInsertSize();
          if (mCurrentVariables.mMaxInsertSize < insertSize) {
            mCurrentVariables.mMaxInsertSize = insertSize;
          }
          if (mCurrentVariables.mMinInsertSize > Math.abs(insertSize)) {
            mCurrentVariables.mMinInsertSize = Math.abs(insertSize);
          }
          if (variables.mInsertSizes != null && (samRec.getAlignmentStart() < samRec.getMateAlignmentStart() || (samRec.getAlignmentStart() == samRec.getMateAlignmentStart() && samRec.getFirstOfPairFlag()))) {
            Integer iscount;
            if (!variables.mInsertSizes.containsKey(insertSize)) {
              iscount = 1;
            } else {
              iscount = variables.mInsertSizes.get(insertSize) + 1;
            }
            variables.mInsertSizes.put(insertSize, iscount);
          }
        }
      }
    } // end of SAM file
    if (currTemplate != null) {
      accumulatePileUp(pileUp, currTemplate.length);
    }
    return countPerRead;
  }

  private void checkAlignmentScore(SAMRecord samRec, SortedMap<Integer, Integer> alignmentScoreMap, int computedScore) {
    final Integer as = samRec.getIntegerAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE);
    if (alignmentScoreMap != null) {
      if (as == null) {
        if (!samRec.getReadUnmappedFlag()) {
          mErr.println("Record has no alignment score: " + samRec.getSAMString().trim());
        }
      } else {
        final Integer ascount;
        if (!alignmentScoreMap.containsKey(as)) {
          ascount = 1;
        } else {
          ascount = alignmentScoreMap.get(as) + 1;
        }
        alignmentScoreMap.put(as, ascount);
      }
    }
    if (mPenaltiesSet && computedScore != Integer.MIN_VALUE && as != null && !as.equals(computedScore)) {
      final String superCigar = samRec.getStringAttribute(SamUtils.CG_SUPER_CIGAR);
      if (superCigar == null) {
        mErr.println("Record's " + SamUtils.ATTRIBUTE_ALIGNMENT_SCORE + ": " + as + " disagreed with computed score: " + computedScore + ", " + samRec.getSAMString().trim());
      }
    }
  }

  private void accumulateAlignmentCounts(SAMRecord samRec, SortedMap<Integer, Integer> map, String attribute) {
    if (map != null) {
      final Integer recordCount = samRec.getIntegerAttribute(attribute);
      if (recordCount != null) {
        Integer mapCount;
        if (!map.containsKey(recordCount)) {
          mapCount = 1;
        } else {
          mapCount = map.get(recordCount) + 1;
        }
        map.put(recordCount, mapCount);
      }
    }
  }

  private void accumulatePileUp(final PileUp pileUp, final int templateLength) {
    if (pileUp != null && mShowConsensus) {
      mCurrentVariables.mConsensus += pileUp.consensus();
      mCurrentVariables.mTotalNt += pileUp.total();
      mCurrentVariables.mTotalLength += templateLength;
    }
  }

  private void printBasicStats(final long totalRecords, final long unmappedRecords, final long formatProblems) {
    mOut.println("Total records: " + totalRecords);
    mOut.println("Unmapped records: " + unmappedRecords);
    if (formatProblems != 0) {
      mOut.println("Records failing samtools validation because of invalid SAM format: " + formatProblems);
    }
  }

  private void printStats(SamStatsVariables variables) {
    if (variables.mTotalLength != 0) {
      mOut.println("Coverage: " + Utils.realFormat(variables.mTotalNt / (double) variables.mTotalLength, 4));
    }
    if (variables.mTotalNt != 0) {
      mOut.println("Consensus: " + Utils.realFormat(variables.mConsensus / (double) variables.mTotalNt, 4));
    }
    if (variables.mTotalMatches != 0) {
      mOut.println("% accuracy compared to reference: " + Utils.realFormat(variables.mTotalMatches / (double) (variables.mTotalMatches + variables.mTotalMismatches) * 100, 2));
    }
    if (!mSingleEnded && variables.mPairOrientations[0] + variables.mPairOrientations[1] + variables.mPairOrientations[2] + variables.mPairOrientations[3] > 0) {
      mOut.println();
      mOut.println("FF: " + variables.mPairOrientations[0]);
      mOut.println("RF: " + variables.mPairOrientations[1]);
      mOut.println("FR: " + variables.mPairOrientations[2]);
      mOut.println("RR: " + variables.mPairOrientations[3]);
    }
    mOut.println();
    if (variables.mMinInsertSize <= variables.mMaxInsertSize) {
      mOut.println("Min/Max insert size: " + variables.mMinInsertSize + "/" + variables.mMaxInsertSize);
    }
    if (variables.mMinReadLength <= variables.mMaxReadLength) {
      mOut.println("Min/Max read length: " + variables.mMinReadLength + "/" + variables.mMaxReadLength);
    }
  }

  private void printDistributions(SamStatsVariables variables, PrintStream out) {
    if (variables.mAlignmentScores != null && variables.mAlignmentScores.size() > 0) {
      printDistribution(variables.mAlignmentScores, "alignment score", "score", "AS", out);
    }
    if (variables.mNH != null && variables.mNH.size() > 0) {
      printDistribution(variables.mNH, "read hits", "hits", "NH", out);
    }
    if (variables.mIH != null && variables.mIH.size() > 0) {
      printDistribution(variables.mIH, "read hits", "hits", "IH", out);
    }
    if (!mSingleEnded && variables.mInsertSizes != null && variables.mInsertSizes.size() > 0) {
      printDistribution(variables.mInsertSizes, "insert size", "size", "IS", out);
    }
  }

  private static void printDistribution(SortedMap<Integer, Integer> values, String titleChunk, String columnChunk, String label, PrintStream out) {
    out.println();
    out.println("Distribution of " + titleChunk + ":");
    out.println("#tag\t" + columnChunk + "\tcount");
    for (final Entry<Integer, Integer> e : values.entrySet()) {
      out.println(label + ":\t" + e.getKey() + "\t" + e.getValue());
    }
  }

  static void printStats2(int[] countPerRead, boolean[] pairedRead, PrintStream out, boolean singleEnded) {
    if (countPerRead != null && ((pairedRead != null && pairedRead.length == countPerRead.length) || singleEnded)) {
      int nonZeroReads = 0;
      int maxCount = 0;
      int uniquelyMapped = 0;
      final TreeMap<Integer, Integer> countsByCount = new TreeMap<>();
      for (int i = 0; i < countPerRead.length; i++) {
        final int count = countPerRead[i];
        final boolean paired = (!(singleEnded || pairedRead == null)) && pairedRead[i];
        if (!countsByCount.containsKey(count)) {
          countsByCount.put(count, 1);
        } else {
          final Integer mapValue = countsByCount.get(count);
          countsByCount.put(count, mapValue + 1);
        }

        if (count > 0) {
          nonZeroReads++;
          if (count > maxCount) {
            maxCount = count;
          }
          if (singleEnded) {
            if (count == 1) {
              uniquelyMapped++;
            }
          } else {
            if (count == 2 && paired) {
              uniquelyMapped++;
            }
          }
        }
      }
      out.println("Reads seen: " + nonZeroReads + " of " + countPerRead.length);
      out.println("Reads mapped at a single location: " + uniquelyMapped);
      out.println("Max records for a read: " + maxCount);
      out.println();
      out.println("Record counts by read:");
      for (final Entry<Integer, Integer> entry : countsByCount.entrySet()) {
        out.println(entry.getKey() + " : " + entry.getValue());
      }
    }
  }

  private boolean matchesGotohCg(byte[] read, byte[] quality, SAMRecord record) {
    mSuperCigarValidator.setTemplateStart(record.getAlignmentStart() - 1);
    try {
      mSuperCigarValidator.setData(record, read, quality);
      mSuperCigarValidator.parse();
    } catch (final IllegalStateException ise) {
      mErr.println(ise.getMessage());
      return false;
    } catch (final BadSuperCigarException bce) {
      mErr.println(bce.getMessage() + ", " + record.getSAMString().trim());
      return false;
    }
    if (!mSuperCigarValidator.isValid()) {
      mErr.println(mSuperCigarValidator.getInvalidReason());
      return false;
    }

    return true;
  }

  /**
   * Check the read as represented in the sam file against the raw read
   * @param read raw read data
   * @param quality raw quality data
   * @param record sam record
   * @param cgData true if CG data
   * @param ignoreFragmentSize true to ignore the CG fragment size
   * @return true if data match
   */
  boolean matchesRawRead(byte[] read, byte[] quality, SAMRecord record, boolean cgData, boolean ignoreFragmentSize) {
    if (read == null) {
      return true;
    } else if (cgData) {
      final String superCigar = record.getStringAttribute(SamUtils.CG_SUPER_CIGAR);
      if (superCigar != null) {
        return matchesGotohCg(read, quality, record);
      }
      return SamValidatorCgHelper.matchesCg(read, quality, record, ignoreFragmentSize);
    }
    final byte[] recordBytes = record.getReadBases();
      if (read.length != recordBytes.length) {
        return false;
      } else if (quality != null && quality.length != record.getBaseQualities().length) {
        return false;
      }
    if (record.getReadNegativeStrandFlag()) {
      for (int i = 0; i < read.length; i++) {

          if (DnaUtils.getBase(DNA.complement(read[i])) != Character.toUpperCase((char) recordBytes[recordBytes.length - i - 1])) {
            return false;
          } else if (quality != null && i < record.getBaseQualities().length && quality[i] != record.getBaseQualities()[recordBytes.length - i - 1]) {
            return false;
          }
      }
    } else {
      for (int i = 0; i < read.length; i++) {
        if (DnaUtils.getBase(read[i]) != Character.toUpperCase((char) recordBytes[i])) {
          return false;
        } else if (quality != null && i < record.getBaseQualities().length && quality[i] != record.getBaseQualities()[i]) {
          return false;
        }
      }
    }
    return true;
  }

  /**
   * Check whether the sam record is valid against the template
   * @param template the template bytes
   * @param samRecord the sam record to check
   * @param pileUp Pile Up object to accumulate consensus statistics to
   * @return -1 if the sam record is invalid, integer min value if score should be ignored, otherwise the computed alignment score.
   */
  public int isAtExpectedRef(final byte[] template, final SAMRecord samRecord, final PileUp pileUp) {
    if (samRecord.getAlignmentStart() < 1 || samRecord.getAlignmentStart() > template.length) {
      mErr.println("Match start position exceeds template limits");
      return -1;
    }
    final byte[] read = samRecord.getReadBases();
    final String cigar = samRecord.getCigarString();

    int n = 0;
    int rPos = 0;
    int tPos = samRecord.getAlignmentStart() - 1;
    int mismatches = 0;

    char prevChar = (char) -1;
    int score = 0;
    boolean ignoreScore = false;

    for (int i = 0; i < cigar.length(); i++) {
      final char c = cigar.charAt(i);
      if (Character.isDigit(c)) {
        n = 10 * n + c - '0';
      } else {
        assert n > 0;
        for (int j = 0; j < n; j++) {
          if (tPos >= template.length && c != SamUtils.CIGAR_SOFT_CLIP) {
            mErr.println("Template length exceeded but read does not indicate soft clipping, " + samRecord.toString());
            return -1;
          }
          if (c == SamUtils.CIGAR_SAME_OR_MISMATCH) { //match OR mismatch
            if (rPos >= read.length) {
              mErr.println("Match went off end of read, " + samRecord.toString());
              return -1;
            }
            final int nt = Character.toUpperCase((char) read[rPos]);
            final int refNt = DnaUtils.getBase(template[tPos]);
            mCurrentVariables.mTotalMatches++;
            if (pileUp != null) {
              pileUp.add((char) nt, tPos);
            }
            if (nt == 'N' || refNt == 'N') { //if either are unknown
              mismatches++;
              mCurrentVariables.mTotalMismatches++;
              score += mUnknownsPenalty;
            } else if (!(nt == refNt)) { // if nt are different
              mismatches++;
              mCurrentVariables.mTotalMismatches++;
              score += mMismatchPenalty;
            }
            tPos++;
            rPos++;
          } else if (c == SamUtils.CIGAR_SAME) { //match
            if (rPos >= read.length) {
              mErr.println("Match went off end of read, " + samRecord.toString());
              return -1;
            }
            final int nt = Character.toUpperCase((char) read[rPos]);
            final int refNt = DnaUtils.getBase(template[tPos]);
            mCurrentVariables.mTotalMatches++;
            if (pileUp != null) {
              pileUp.add((char) nt, tPos);
            }
            if (nt == 'N' || refNt == 'N') {
              score += mUnknownsPenalty;
            } else if (nt != refNt) { // if nt are different or, both are N, or either are N, mismatch.
              mErr.println("Expected match " + (char) refNt + " was " + (char) nt + ", rpos=" + rPos + ", " + samRecord.toString());
              return -1;
            }
            tPos++;
            rPos++;
          } else if (c == SamUtils.CIGAR_MISMATCH) { //mismatch
            if (rPos >= read.length) {
              mErr.println("Match went off end of read, " + samRecord.toString());
              return -1;
            }
            final int nt = Character.toUpperCase((char) read[rPos]);
            final int refNt = DnaUtils.getBase(template[tPos]);
            mCurrentVariables.mTotalMatches++;
            if (pileUp != null) {
              pileUp.add((char) nt, tPos);
            }

            if (nt == 'N' || refNt == 'N') {
              score += mUnknownsPenalty;
            } else if (nt != refNt) {
              score += mMismatchPenalty;
            } else {
              mErr.println("Expected mismatch " + samRecord.toString());
              return -1;
            }
            mCurrentVariables.mTotalMismatches++;
            mismatches++;
            tPos++;
            rPos++;
          } else if (c == SamUtils.CIGAR_DELETION_FROM_REF) {
            tPos++;
            mismatches++;
            if (prevChar != SamUtils.CIGAR_DELETION_FROM_REF) {
              score += mGapOpenPenalty;
            }
            score += mGapExtendPenalty;
          } else if (c == SamUtils.CIGAR_GAP_IN_READ) { // skip used in CG reads
            tPos++;
          } else if (c == SamUtils.CIGAR_INSERTION_INTO_REF) {
            mismatches++;
            rPos++;
            if (prevChar != SamUtils.CIGAR_INSERTION_INTO_REF) {
              score += mGapOpenPenalty;
            }
            score += mGapExtendPenalty;
          } else if (c == SamUtils.CIGAR_SOFT_CLIP) { // soft-clipping bases in read ignored for position
            rPos++; // NM field does not count soft clip, hence don't increment mismatches
            ignoreScore = true;
          } else {
            assert false;
          }
          prevChar = c;
        }
        n = 0;
      }
    }
    boolean ok = true;
    if (samRecord.getIntegerAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES) != null) {
      final Integer samnm = samRecord.getIntegerAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES);
      if (samnm == null) {
        mErr.println("SAM record did not contain " + SamUtils.ATTRIBUTE_NUM_MISMATCHES + " attribute. " + samRecord.toString());
        return -1;
      }
      ok = samnm == mismatches;
    }
    if (!ok) {
      dumpFaulty(template, samRecord, mErr); // Attempt to print the faulty alignment
      mErr.println("Observed mismatches: " + mismatches + " Claimed mismatches: " + samRecord.getIntegerAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES));
    }
    return ignoreScore ? Integer.MIN_VALUE : ok ? score : -1;
  }

  private static char readNt(final byte[] read, final int p) {
    return Character.toUpperCase(p < read.length ? (char) read[p] : '#');
  }

  private static char refNt(final byte[] seq, final int p) {
    return DnaUtils.getBase(seq[p]);
  }

  private static void dumpFaulty(final byte[] sequence, final SAMRecord samRecord, final PrintStream out) {
    final byte[] read = samRecord.getReadBases();
    final String cigar = samRecord.getCigarString();
    int n = 0;
    int rPos = 0;
    int tPos = samRecord.getAlignmentStart() - 1;
    final StringBuilder sbRead = new StringBuilder();
    final StringBuilder sbMatches = new StringBuilder();
    final StringBuilder sbTemplate = new StringBuilder();
    for (int i = 0; i < cigar.length(); i++) {
      final char c = cigar.charAt(i);
      if (Character.isDigit(c)) {
        n = 10 * n + c - '0';
      } else {
        for (int j = 0; j < n && tPos < sequence.length; j++) {
          switch (c) {
            case SamUtils.CIGAR_SAME_OR_MISMATCH: //match OR mismatch
            case SamUtils.CIGAR_SAME: //match
            case SamUtils.CIGAR_MISMATCH: //mismatch
              final char nt = readNt(read, rPos++);
              final char refNt = refNt(sequence, tPos++);
              sbRead.append(nt);
              sbMatches.append(nt != refNt && nt != 'N' && refNt != 'N' ? ' ' : '|');
              sbTemplate.append(refNt);
              break;
            case SamUtils.CIGAR_DELETION_FROM_REF:
              sbRead.append('-');
              sbMatches.append(' ');
              sbTemplate.append(refNt(sequence, tPos++));
              break;
            case SamUtils.CIGAR_GAP_IN_READ: // skip used in CG reads
              sbRead.append('.');
              sbMatches.append(' ');
              sbTemplate.append(refNt(sequence, tPos++));
              break;
            case SamUtils.CIGAR_INSERTION_INTO_REF:
              sbTemplate.append('-');
              sbRead.append(readNt(read, rPos++));
              sbMatches.append(' ');
              break;
            case SamUtils.CIGAR_SOFT_CLIP:
              // soft-clipping bases in read ignored for position
              sbRead.append(readNt(read, rPos++));
              sbMatches.append(' ');
              sbTemplate.append(' ');
              break;
            default:
              assert false;
          }
        }
        n = 0;
      }
    }
    if (sbRead.length() < 120) {
      out.println("Claimed alignment is incorrect:");
      out.println(sbTemplate.toString() + " <-- template");
      out.println(sbMatches.toString());
      out.println(sbRead.toString() + " <-- read");
    }
  }

  private byte[] read(final int readId, final boolean first) {
    return ReadHelper.getRead(first ? mLeftReader : mRightReader, readId);
  }

  private byte[] qual(final int readId, final boolean first) {
    return ReadHelper.getQual(first ? mLeftReader : mRightReader, readId);
  }

  static final class SamStatsVariables {
    int mMaxInsertSize = 0;
    int mMinInsertSize = Integer.MAX_VALUE;
    int mMaxReadLength = 0;
    int mMinReadLength = Integer.MAX_VALUE;
    long mTotalRecords = 0;
    long mTotalMatches = 0;
    long mTotalMismatches = 0;
    long mUnmappedRecords = 0;
    final int[] mPairOrientations = new int[4];  //0=FR, 1=RF, 2=FF, 3=RR
    long mConsensus = 0;
    long mTotalNt = 0;
    long mTotalLength = 0;
    final SortedMap<Integer, Integer> mAlignmentScores = new TreeMap<>();
    final SortedMap<Integer, Integer> mInsertSizes = new TreeMap<>();
    final SortedMap<Integer, Integer> mNH = new TreeMap<>();
    final SortedMap<Integer, Integer> mIH = new TreeMap<>();

    void addToTotal(SamStatsVariables individual) {
      if (mMaxInsertSize < individual.mMaxInsertSize) {
        mMaxInsertSize = individual.mMaxInsertSize;
      }
      if (mMinInsertSize > individual.mMinInsertSize) {
        mMinInsertSize = individual.mMinInsertSize;
      }
      if (mMaxReadLength < individual.mMaxReadLength) {
        mMaxReadLength = individual.mMaxReadLength;
      }
      if (mMinReadLength > individual.mMinReadLength) {
        mMinReadLength = individual.mMinReadLength;
      }
      mTotalRecords += individual.mTotalRecords;
      mTotalMatches += individual.mTotalMatches;
      mTotalMismatches += individual.mTotalMismatches;
      mUnmappedRecords += individual.mUnmappedRecords;
      for (int i = 0; i < mPairOrientations.length; i++) {
        mPairOrientations[i] += individual.mPairOrientations[i];
      }
      mConsensus += individual.mConsensus;
      mTotalNt += individual.mTotalNt;
      mTotalLength += individual.mTotalLength;
      mergeMaps(mAlignmentScores, individual.mAlignmentScores);
      mergeMaps(mInsertSizes, individual.mInsertSizes);
      mergeMaps(mNH, individual.mNH);
      mergeMaps(mIH, individual.mIH);
    }

    private void mergeMaps(SortedMap<Integer, Integer> target, SortedMap<Integer, Integer> individual) {
      for (final Entry<Integer, Integer> pair : individual.entrySet()) {
        if (target.containsKey(pair.getKey())) {
          target.put(pair.getKey(), target.get(pair.getKey()) + pair.getValue());
        } else {
          target.put(pair.getKey(), pair.getValue());
        }
      }
    }
  }
}
