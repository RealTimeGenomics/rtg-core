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
package com.rtg.simulation.cnv;

import java.io.IOException;
import java.io.OutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import com.rtg.mode.DNA;
import com.rtg.reader.SdfWriter;
import com.rtg.reader.SequencesReader;
import com.rtg.util.Environment;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.MathUtils;
import com.rtg.util.PortableRandom;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

/**
 */
public class CnvSimulator {
  //TODO: large inserts?
  //TODO: implement some/more priors

  private static final String CNVSIM_VERSION = "v0.0";
  private static final String TB = "\t";
  private final SequencesReader mInput;
  private final SdfWriter mOutput;
  private final SdfWriter mTwin;
  private final OutputStream mCnvOutput;
  private final PortableRandom mRandom;
  private final CnvPriorParams mPriors;
  private final double mPercentCnved;
  protected int mSetCountCnved;

  protected ArrayList<RegionDummy> mRegionDummyList = new ArrayList<>();

  long mTotalSequenceLength;
  int[] mSequenceLengths;
  String[] mSequenceNames;
  long mNumberSequences;
  int[] mCountsPerSeq;
  long[] mLengthCnvedPerSeq;
  int mTotalCnvedCount;
  long mTotalCnvedLength;
  int mTries = 0;
  private int mMinCnvLength;

  protected long mSetCnvedSequenceLength;
  protected long[] mBreakPoints;
  protected List<List<CnvRegion>> mRegionSequences;

  private final double[] mPowerLengthThresholds;
  protected int[] mLengthHistogramCnved = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  protected int[] mLengthHistogramUncnved = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  /**
   * Constructor
   * @param input SDF reader
   * @param output SDF writer for genome
   * @param twin SDF writer for twin genome
   * @param mappingOutput file with information about generated CNVs
   * @param random random number generator
   * @param priors properties file for probabilities and thresholds
   * @param percentCnved percent of the genome that should be part of a CNV region
   * @param countCnved count of CNV regions
   */
  public CnvSimulator(final SequencesReader input, final SdfWriter output, final SdfWriter twin, final OutputStream mappingOutput,
      final PortableRandom random, final CnvPriorParams priors, final double percentCnved, final int countCnved) {

    mInput = input;
    mOutput = output;
    mTwin = twin;
    mCnvOutput = mappingOutput;
    mRandom = random;
    mPriors = priors;
    if (mPriors != null) {
      mMinCnvLength = 0;
    }
    mPercentCnved = percentCnved;
    mSetCountCnved = countCnved;
    mPowerLengthThresholds = mPriors != null ? mPriors.powerLengthThresholds() : null;
  }

  protected static class FixedRegion {
    final String mSequenceName;
    final long mStartPosNullBased;
    final long mEndPosNullBased;
    final int mLength;
    final int mCopyNumber;

    FixedRegion(final String sequenceName, final long startPos, final long endPos, final int copyNumber) {
      if (startPos < 1 || endPos < startPos) {
        throw new InvalidParamsException(ErrorType.INVALID_PARAMETER, "CNV invalid start end positions, start: " + startPos + " end: " + endPos);
      }
      mSequenceName = sequenceName;
      mStartPosNullBased = startPos - 1;
      mEndPosNullBased = endPos - 1;
      mLength = (int) mEndPosNullBased - (int) mStartPosNullBased + 1;
      mCopyNumber = copyNumber;
    }
  }

  protected static FixedRegionComparator getFixedRegionComparator() {
    return new FixedRegionComparator();
  }

  protected static class FixedRegionComparator implements Comparator<FixedRegion>, Serializable {
    @Override
    public int compare(final FixedRegion o1, final FixedRegion o2) {
      return (Integer.valueOf((int) o1.mStartPosNullBased)).compareTo((int) o2.mStartPosNullBased);
    }
  }

  protected static class RegionDummy {
    int mSequenceId;
    int mRegionPos;

    RegionDummy(final int seq, final int pos) {
      mSequenceId = seq;
      mRegionPos = pos;
    }
  }

  protected void generate() throws IOException {

    assert mPriors != null;
    initialize();
    mMinCnvLength = -1;
    mBreakPoints = taskaGetBreakPoints(mTotalSequenceLength, mInput.maxLength(), -1, -1, -1);
    middleTasks();

    endTasks();
  }

  protected void initialize() throws IOException {
    mRegionSequences = new ArrayList<>();
    mTotalSequenceLength = mInput.totalLength();
    mNumberSequences = mInput.numberSequences();
    mSequenceLengths = mInput.sequenceLengths(0, mNumberSequences);
    mSequenceNames = new String[(int) mNumberSequences];
    mCountsPerSeq = new int[(int) mNumberSequences];
    mLengthCnvedPerSeq = new long[(int) mNumberSequences];
    mTotalCnvedCount = 0;
    mTotalCnvedLength = 0;
    for (int i = 0; i < mInput.numberSequences(); i++) {
      mSequenceNames[i] = mInput.name(i);
    }
  }

  protected void generate(final String[] lengthAndVariance) throws IOException, InvalidParamsException {
    initialize();
    try {
      final int avgLength = Integer.parseInt(lengthAndVariance[0]);
      final int lengthDiffBound = Integer.parseInt(lengthAndVariance[1]);
      if (avgLength < 1 || lengthDiffBound < 0) {
        throw new InvalidParamsException(ErrorType.INVALID_PARAMETER, "CNV invalid length or bound for variant: "
            + avgLength + ", " + lengthDiffBound);
      }
      final int rangeLength = lengthDiffBound * 2 + 1;
      mMinCnvLength = avgLength - lengthDiffBound;
      mBreakPoints = taskaGetBreakPoints(mTotalSequenceLength, mInput.maxLength(), lengthDiffBound, avgLength, rangeLength);
    } catch (final NumberFormatException e) {
      throw new InvalidParamsException(ErrorType.INVALID_PARAMETER, "CNV length should contain two integers ");
    }
    middleTasks();
    endTasks();
  }

  private void middleTasks() {
    // generate CNV region objects
    taskbGenerateCnvRegionsObjects();
    // set how many are selected
    taskcSetLengthsAndCount(mMinCnvLength);
    // select for copying
    // select per count of per percent
    taskdSelectCnvs();
    taskePrepareCopying();
  }

  private void endTasks() throws IOException {
    // output to mut. genome and mut. twin
    taskfOutputToGenomeSdfs();
    // accumulate CNV regions
    taskgAccumulateNonCnvRegions();
    // output cnv info file
    outputBreakPointsToInfoFile();
  }

  protected void generate(final FixedRegion[] fixedRegions) throws IOException, InvalidParamsException {

    assert mPriors == null;
    initialize();
    // no breakpoints needed.
    mBreakPoints = null;
    final FixedRegion[][] fixedRegionLists = generatedSortedRegionLists(fixedRegions, (int) mNumberSequences);
    // generate the CNV region and all CNV region objects needed
    taskbGenerateCnvRegionsObjects(fixedRegionLists);
    endTasks();
  }

  private FixedRegion[][] generatedSortedRegionLists(final FixedRegion[] fixedRegions, final int numberSequences) {
    @SuppressWarnings("unchecked")
    final ArrayList<FixedRegion>[] list = (ArrayList<FixedRegion>[]) new ArrayList<?>[numberSequences];
    for (int i = 0; i < list.length; i++) {
      list[i] = new ArrayList<>();
    }
    for (final FixedRegion region : fixedRegions) {
      boolean found = false;
      for (int i = 0; i < list.length; i++) {
        if (mSequenceNames[i].equals(region.mSequenceName)) {
          found = true;
          list[i].add(region);
          break;
        }
      }
      if (!found) {
        throw new InvalidParamsException(ErrorType.INVALID_PARAMETER, "Wrong sequence name.");
      }
    }
    final FixedRegion[][] sorted = new FixedRegion[numberSequences][];
    for (int i = 0; i < numberSequences; i++) {
      sorted[i] = new FixedRegion[list[i].size()];
      for (int j = 0; j < list[i].size(); j++) {
        sorted[i][j] = list[i].get(j);
      }
      Arrays.sort(sorted[i], getFixedRegionComparator());
    }
    return sorted;
  }

  protected void taskbGenerateCnvRegionsObjects() {
    long currentSequenceStart = 0;
    int currentBPIndex = 0;
    long lastBreakPoint = 0;
    for (int i = 0; i < mNumberSequences; i++) {
      final ArrayList<CnvRegion> list = new ArrayList<>();
      while (currentBPIndex < mBreakPoints.length
          && mBreakPoints[currentBPIndex] - currentSequenceStart < mSequenceLengths[i]) {

        final CnvRegion region = new CnvRegion(i, (int) (lastBreakPoint - currentSequenceStart),
            (int) (mBreakPoints[currentBPIndex] - lastBreakPoint), mPriors);
        // add region
        mRegionDummyList.add(new RegionDummy(i, list.size()));
        list.add(region);

        lastBreakPoint = mBreakPoints[currentBPIndex];
        currentBPIndex++;
      }
      final int start = (int) (lastBreakPoint - currentSequenceStart);
      // add last region but not end region
      mRegionDummyList.add(new RegionDummy(i, list.size()));
      list.add(new CnvRegion(i, start, mSequenceLengths[i] - start, mPriors));
      list.add(new CnvRegion(i, mSequenceLengths[i], 0)); // endRegion
      currentSequenceStart += mSequenceLengths[i];
      lastBreakPoint = currentSequenceStart;
      mRegionSequences.add(list);
    }
  }

  protected void taskbGenerateCnvRegionsObjects(final FixedRegion[][] fixedRegions) {
    for (int i = 0; i < mInput.numberSequences(); i++) {
      final ArrayList<CnvRegion> list = new ArrayList<>();
      //      mInput.seek(i);
      //      final String seqName = mInput.currentName();
      //      mSequenceNames[i] = seqName;
      if (fixedRegions[i].length > 0) {
        // first gap
        if (fixedRegions[i][0].mStartPosNullBased > 0) {
          final CnvRegion startRegion = new CnvRegion(i, 0, (int) fixedRegions[i][0].mStartPosNullBased);
          list.add(startRegion);
        }

        for (int j = 0; j < fixedRegions[i].length; j++) {
          final FixedRegion currentFixed = fixedRegions[i][j];
          if ((int) currentFixed.mStartPosNullBased >= mSequenceLengths[i]) {
            throw new InvalidParamsException(ErrorType.INVALID_PARAMETER, "CNV start too large " + (currentFixed.mStartPosNullBased + 1)
                + " sequence length: " + mSequenceLengths[i]);
          }

          // CNV
          final int posAfter = j + 1 >= fixedRegions[i].length ? mSequenceLengths[i]
                                                                                        : (int) fixedRegions[i][j + 1].mStartPosNullBased;
          if ((int) currentFixed.mEndPosNullBased >= posAfter) {
            throw new InvalidParamsException(ErrorType.INVALID_PARAMETER, "CNV end " + (posAfter + 1)
                + " overlaps or is longer than sequence. Sequence name/length " + mSequenceNames[i]
                                                                                                 + "/" + mSequenceLengths[i]);
          }

          final CnvRegion regionCnved = new CnvRegion(i, (int) currentFixed.mStartPosNullBased,
              currentFixed.mLength, currentFixed.mCopyNumber);
          addCountStatistics(regionCnved);
          list.add(regionCnved);

          CnvRegion goalRegion = regionCnved;

          if (currentFixed.mEndPosNullBased < posAfter - 1) {
            // gap
            final CnvRegion regionGap = new CnvRegion(i, (int) (currentFixed.mEndPosNullBased + 1),
                (int) (posAfter - currentFixed.mEndPosNullBased - 1));
            list.add(regionGap);
            // copy to end of gap
            goalRegion = regionGap;
          }

          // prepare copies
          if (currentFixed.mCopyNumber > 2) {
            for (int k = 2; k < currentFixed.mCopyNumber; k++) {
              goalRegion.addCopy(regionCnved);
            }
          }
        }
      } else {
        final CnvRegion onlyRegion = new CnvRegion(i, 0, mSequenceLengths[i]);
        list.add(onlyRegion);
      }
      // endRegion
      list.add(new CnvRegion(i, mSequenceLengths[i], 0));
      mRegionSequences.add(list);
    }
  }

  private void addCountStatistics(final CnvRegion regionCnved) {
    mTotalCnvedCount++;
    mTotalCnvedLength += regionCnved.mLength;
    mCountsPerSeq[regionCnved.mSequenceId]++;
    mLengthCnvedPerSeq[regionCnved.mSequenceId] += regionCnved.mLength;
    addToHistogram(regionCnved, mLengthHistogramCnved);
  }

  /**
   * select all numbers of copies and positions where to copy to
   */
  protected void taskdSelectCnvs() {
    CnvRegion pickedRegion;
    do {
      mTries = 0;
      do {
        mTries++;
        final RegionDummy randomPick = mRegionDummyList.get(mRandom.nextInt(mRegionDummyList.size()));
        pickedRegion = mRegionSequences.get(randomPick.mSequenceId).get(randomPick.mRegionPos);
        if (mTries > 1000) {
          break;
        }
      } while (pickedRegion.isUnpickable(mSequenceLengths, true, 0));
      if (mTries > 1000) {

        break;
      }
      // select region as CNV region
      pickedRegion.initializeAsCnved(mRandom);

      addCountStatistics(pickedRegion);
    } while(mTotalCnvedCount < mSetCountCnved && mTotalCnvedLength < mSetCnvedSequenceLength);

  }

  protected void taskePrepareCopying() {
    // prepare copying
    for (int i = 0; i < mRegionSequences.size(); i++) {
      for (int j = 0; j < mRegionSequences.get(i).size(); j++) {
        final CnvRegion currentRegion = mRegionSequences.get(i).get(j);
        for (int k = 0; k < currentRegion.mNumCopies; k++) {
          final int destinationId;
          if (mRandom.nextBoolean()) {
            // copied to same sequence
            destinationId = i;
          } else {
            //copied to another sequence
            destinationId = mRandom.nextInt(mRegionSequences.size());
          }
          final List<CnvRegion> destinationSequence = mRegionSequences.get(destinationId);
          final CnvRegion destinationRegion = destinationSequence.get(mRandom.nextInt(destinationSequence.size()));
          boolean isHetero = false;
          if (currentRegion.mNumCopies - k > 1) {
            // enough for 2 copies
            if (mRandom.nextBoolean()) {
              k++;
              isHetero = true;
            }
          }
          destinationRegion.addCopy(currentRegion, mRandom, isHetero);
        }
      }
    }
  }

  protected void taskcSetLengthsAndCount(int minLength) {
    if (mSetCountCnved < Integer.MAX_VALUE) {
      if (mRegionDummyList == null) {
        throw new NoTalkbackSlimException("Cannot find enough regions for " + mSetCountCnved + " CNV regions; 0 possible.");

      } else if (mRegionDummyList.size() <= mSetCountCnved) {
        int possNum = 0;
        for (final RegionDummy regDummy : mRegionDummyList) {
          final CnvRegion region = mRegionSequences.get(regDummy.mSequenceId).get(regDummy.mRegionPos);
          if (!region.isUnpickable(mSequenceLengths, false, minLength)) {
            possNum++;
          }
        }
        throw new NoTalkbackSlimException("Cannot find enough regions for " + mSetCountCnved + " CNV regions; "
            + possNum + " possible.");
      }
      mSetCnvedSequenceLength = Long.MAX_VALUE;
    } else {
      mSetCnvedSequenceLength = MathUtils.round((double) mTotalSequenceLength * (mPercentCnved / 100.0));
    }
  }

  /**
   * Output bases to both genomes, copies region according to number of copies selected
   */
  private void taskfOutputToGenomeSdfs() throws IOException {
    for (int i = 0; i < mNumberSequences; i++) {
      final byte[] genomeSeq = new byte[mInput.length(i)];
      mInput.read(i, genomeSeq);
      //System.err.println("seek "+ i);
      mOutput.startSequence(mSequenceNames[i]);
      mTwin.startSequence(mSequenceNames[i]);
      for (int j = 0; j < mRegionSequences.get(i).size(); j++) {
        final CnvRegion currentRegion = mRegionSequences.get(i).get(j);

        // copy other regions before region
        for (int k = 0; k < currentRegion.mCopies.size(); k++) {
          final byte[] copyDnas = currentRegion.getCopy(k, mInput);
          final boolean bothStrands = currentRegion.getBothStrands(k);
          if (mRandom.nextBoolean()) {
            // reverse
            DNA.reverseComplementInPlace(copyDnas, 0, copyDnas.length);
          }
          if (bothStrands) {
            writeDnaArray(mOutput, copyDnas);
            writeDnaArray(mTwin, copyDnas);
            //currentRegion.mCopies.get(k).mNumCopies++;
          } else {
            if (mRandom.nextBoolean()) {
              writeDnaArray(mOutput, copyDnas);
            } else {
              writeDnaArray(mTwin, copyDnas);
            }
          }
        }

        // write current region
        if (currentRegion.getCN() > 0) {
          final byte[] currentDnas = getDnaArray(genomeSeq, currentRegion.mStart, currentRegion.mLength);
          if (!currentRegion.mOutputDelete) {
            writeDnaArray(mOutput, currentDnas);
          }
          if (!currentRegion.mTwinDelete) {
            writeDnaArray(mTwin, currentDnas);
          }
        }

      }
      mOutput.endSequence();
      mTwin.endSequence();
    }
  }


  /**
   * accumulate non CNV regions together, must happen after writing to genome SDFs
   */
  protected void taskgAccumulateNonCnvRegions() {
    for (int i = 0; i < mInput.numberSequences(); i++) {
      // don't swallow end break point -> size - 1
      //System.err.println("numBP " + mRegionSequences.get(i).size());
      int sizeWithoutEnd = mRegionSequences.get(i).size() - 1;
      for (int j = 0; j < sizeWithoutEnd; j++) {
        final CnvRegion currentRegion = mRegionSequences.get(i).get(j);
        if (currentRegion.mIsCnved) {
          continue;
        } else {
          final int nextAfter = j + 1;
          if (nextAfter < sizeWithoutEnd) {
            CnvRegion testedRegion = mRegionSequences.get(i).get(nextAfter);
            while (!testedRegion.mIsCnved && nextAfter < sizeWithoutEnd) {
              currentRegion.addRegion(testedRegion);
              mRegionSequences.get(i).remove(nextAfter);
              // size changes
              sizeWithoutEnd = mRegionSequences.get(i).size() - 1;
              if (nextAfter < sizeWithoutEnd) {
                testedRegion = mRegionSequences.get(i).get(nextAfter);
              }
            }
          }
          addToHistogram(currentRegion, mLengthHistogramUncnved);
        }
      }
    }
  }

  private void outputBreakPointsToInfoFile() throws IOException {
    try (OutputStream cnvOutput = mCnvOutput) {
      cnvOutput.write(headerLines().getBytes());
      for (int i = 0; i < mRegionSequences.size(); i++) {
        for (int j = 0; j < mRegionSequences.get(i).size(); j++) {
          final CnvRegion currentRegion = mRegionSequences.get(i).get(j);
          //System.err.println(currentRegion.toString());
          // write info file about generated CNV
          final String seqName = mInput.name(i);
          cnvOutput.write(currentRegion.toBytes(i, seqName));
        }
      }
    }
  }

  /**
   * Prints the histograms for the region lengths
   * @param foStream out put stream
   * @throws IOException when writing goes wrong
   */
  public void outputHistograms(final OutputStream foStream) throws IOException {
    foStream.write(StringUtils.LS.getBytes());

    if (mTries > 1000) {
      foStream.write(("Gave up after " + mTries + " attempts.").getBytes());
    }
    foStream.write(StringUtils.LS.getBytes());

    //System.err.println("nts " + mTotalCnvedLength);
    foStream.write(("Total length" + TB + mTotalSequenceLength + TB
        + "CNV-count" + TB + mTotalCnvedCount + TB + "CNV-percent" + TB
        +  percentString(mTotalCnvedLength, mTotalSequenceLength)).getBytes());

    foStream.write(StringUtils.LS.getBytes());

    for (int i = 0; i < mNumberSequences; i++) {
      //System.err.println("nts[" + i + "]: " + mLengthPerSeq[i]);
      foStream.write(("Sequence" + TB + i + TB + "length" + TB + mSequenceLengths[i] + TB
          + "CNV-count" + TB + mCountsPerSeq[i] + TB + "CNV-percent" + TB
          + percentString(mLengthCnvedPerSeq[i], mSequenceLengths[i]) + " name: " + mSequenceNames[i]).getBytes());
    }
    //    printStream.println("");
    //    printStream.println("Lengths INBETWEEN CNV regions");
    //    outputHistogram(out, mLengthHistogramUncnved);
    foStream.write(StringUtils.LS.getBytes());
    foStream.write("Lengths of CNV regions".getBytes());
    foStream.write(StringUtils.LS.getBytes());
    outputHistogram(foStream, mLengthHistogramCnved);

    //System.err.println("nts " + mTotalCnvedLength);
    foStream.write(("Total length" + TB + mTotalSequenceLength + StringUtils.LS
        + "CNV-count" + TB + mTotalCnvedCount + StringUtils.LS + "CNV-percent" + TB
        +  percentString(mTotalCnvedLength, mTotalSequenceLength) + StringUtils.LS).getBytes());
  }

  protected static String percentString(final long seqLengthCnved, final long seqLength) {
    return Utils.realFormat(((double) seqLengthCnved / (double) seqLength) * 100, 2) + "%";
  }

  static final String[] HISTOGRAM_TEXT = {
    "[0               - 0] : ",
    "[1              - 10] : ",
    "[11            - 100] : ",
    "[101         - 1,000] : ",
    "[1,001      - 10,000] : ",
    "[1,0001    - 100,000] : ",
    "[..      - 1,000,000] : ",
    "[..     - 10,000,000] : ",
    "[..    - 100,000,000] : ",
  "[..  - 1,000,000,000] : "};

  private void outputHistogram(final OutputStream out, final int[] histo) throws IOException {
    out.write(("range of nt-length    : count" + StringUtils.LS).getBytes());
    for (int i = 0; i < histo.length; i++) {
      out.write((HISTOGRAM_TEXT[i] + histo[i] + StringUtils.LS).getBytes());
    }
  }

  protected String regionSequencesToString() {
    final StringBuilder sb2 = new StringBuilder();
    for (int i = 0; i < mRegionSequences.size(); i++) {
      for (int j = 0; j < mRegionSequences.get(i).size(); j++) {
        final CnvRegion r = mRegionSequences.get(i).get(j);
        sb2.append(r.toString());
      }
    }
    return sb2.toString();
  }

  /**
   * Choose from accumulated distribution
   * @param rand a double chosen between 0.0 and 1.0
   * @param thresholds accumulated distribution values
   * @return a chosen index
   */
  static int chooseFromAccumDistribution(final double rand, final double[] thresholds) {
    int index = 0;
    while (thresholds[index] < rand) {
      index++;
    }
    return index;
  }

  protected void addToHistogram(final CnvRegion currentRegion, final int[] histo) {
    histo[getMagnitudeIndex(currentRegion.mLength, histo.length - 1)]++;
  }

  /**
   * @param value value to find the magnitude of
   * @param maxIndex upper bound on the index
   * @return index according to the magnitude of the value
   */
  static int getMagnitudeIndex(final int value, final int maxIndex) {
    return (value > 0) ? Math.min(maxIndex, 1 + (int) Math.log10(value)) : 0;
  }

  /** Header line for CNV bed file. */
  public static final String CNV_HEADER = "#Seq" + TB + "start" + TB + "end" + TB + "label" + TB + "cn" + TB + "bp-cn" + TB + "error" + StringUtils.LS;

  protected static String headerLines() {
    final StringBuilder line = new StringBuilder();
    line.append("#Version ").append(Environment.getVersion()).append(", cnvsim ").append(CNVSIM_VERSION).append(StringUtils.LS);
    if (CommandLine.getCommandLine() != null) {
      line.append("#CL" + TB).append(CommandLine.getCommandLine()).append(StringUtils.LS);
    } else {
      line.append("#CL" + TB + "Internal").append(StringUtils.LS);
    }
    line.append(CNV_HEADER);
    return line.toString();
  }

  static byte[] getDnaArray(final byte[] seq, final int pos, final int length) {
    //System.err.println("sequence length=" + seq.length + " pos=" + pos + " length=" + length);
    assert pos + length <= seq.length;
    final byte[] arr = new byte[length];
    for (int i = pos, j = 0; i < pos + length; i++, j++) {
      arr[j] = seq[i];
    }
    return arr;
  }

  private void writeDnaArray(final SdfWriter out, final byte[] dnas) throws IOException {
    if (dnas != null) {
      out.write(dnas, null, dnas.length);
    }
  }

  protected long[] taskaGetBreakPoints(final long totalSequenceLength, final long maxLength, int avgLength, int lengthDiffBound, int rangeLength) {
    final ArrayList<Long> bps = new ArrayList<>();
    long currentLength = 0;
    while (currentLength < totalSequenceLength) {
      long length = Long.MAX_VALUE;
      // find one that is not too long to fit anywhere
      // if the distribution contains ranges larger than maximum sequence size, these ranges are ignored
      while (length > maxLength) {
        if (avgLength > 0) {
          // length is given per parameter
          length = lengthDiffBound == 0 ? avgLength : mMinCnvLength + mRandom.nextInt(rangeLength);
        } else {
          final int power = chooseFromAccumDistribution(mRandom.nextDouble(), mPowerLengthThresholds);
          length = randomPriorsLength(mRandom.nextDouble(), power);
        }
      }
      currentLength += length;
      if (currentLength < totalSequenceLength) {
        bps.add(bps.size() == 0 ? 0 : mRandom.nextInt(bps.size()), length);
      }
    }
    final long[] results = new long[bps.size()];
    long position = 0;
    for (int i = 0; i < results.length; i++) {
      position += bps.get(i);
      results[i] = position;
    }
    return results;
  }


  /**
   * @param rand value between 0 and 1
   * @param power of which magnitude the length is
   * @return new length to use for CNV
   */
  protected long randomPriorsLength(final double rand, final int power) {
    return (long) (rand * Math.pow(10, (double) power + 1) + 1);
  }

  @Override
  public String toString() {
    final StringBuilder line = new StringBuilder();
    line.append("CnvSimulator").append(StringUtils.LS);
    line.append("No. breakpoints ").append(mBreakPoints.length).append(StringUtils.LS);
    for (final long bp : mBreakPoints) {
      line.append(bp).append(StringUtils.LS);
    }
    for (int i = 0; i < mRegionSequences.size(); i++) {
      line.append("Sequence ").append(i).append(" No. regions ").append(mRegionSequences.get(i).size()).append(StringUtils.LS);
    }
    if (mPriors != null) {
      line.append("priors set");
    }
    return line.toString();

  }
}
