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
package com.rtg.variant.sv;

import static com.rtg.util.StringUtils.TAB;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import com.rtg.launcher.NoStatistics;
import com.rtg.sam.ReadGroupUtils;
import com.rtg.sam.SamIteratorTask;
import com.rtg.util.Environment;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.util.VariantUtils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;

/**
 */
public class SvToolTask extends SamIteratorTask<SvToolParams, NoStatistics> {

  private static final String SV_OUTPUT_VERSION = "v0.1";

  private static final int MAX_WARNINGS = 5;

  private final Map<String, String> mReadGroupLabels;
  private final Map<String, ReadGroupState> mReadGroupStates;

  private Signal[] mSimpleSignals;
  private Signal[] mPosteriorSignals;

  private OutputStream mSimpleOut;
  private OutputStream mBayesianOut;
  private int mWarnCount = 0;

  private int mGlobalRadius = 0;

  /** The actual bytes of the current template. */
  private byte[] mTemplate;
  private final SamCounts mTemplateNs = new CumulativeSamCounts(0, null);
  private int mTemplatePos;
  private int mTemplateMin;

  private int mPreviousSam;

  private final int mDefaultStepSize;
  private final int mZoomedStepSize;
  private int mCurrentStepSize;
  private int mLastChange = -1;
  private int mLastHypothesis;
  private int mLastWrite = 0;
  private final Corrections mCorrections;

  private SvInterestingRegionExtractor mInterestingRegionExtractor;

  protected SvToolTask(SvToolParams params, OutputStream defaultOutput) throws IOException {
    super(params, defaultOutput, new NoStatistics(), params.filterParams());
    if (mGenomeSequences == null) {
      throw new NoTalkbackSlimException(ErrorType.READING_ERROR, mParams.genome().toString());
    }

    mReadGroupLabels = params.readGroupLabels();
    mReadGroupStates = new HashMap<>();


    mZoomedStepSize = mParams.fineStepSize();
    mDefaultStepSize = mParams.stepSize();

    mCurrentStepSize = mParams.stepSize();
    final File corr = mParams.correctionsFile();
    mCorrections = corr == null ? null : new Corrections(corr);
  }

  private String getReadGroupLabel(String rgId) {
    if (mReadGroupLabels != null && mReadGroupLabels.containsKey(rgId)) {
      return mReadGroupLabels.get(rgId);
    }
    return rgId;
  }

  private String getReadGroupLabel(SAMRecord rec) {
    return getReadGroupLabel(ReadGroupUtils.getReadGroup(rec));
  }

  @Override
  protected void init(SAMFileHeader header) throws IOException {
    final Map<String, ReadGroupStats> rgStats = mParams.readGroupStatistics();

    final Map<String, MachineType> rgMachineTypes = new HashMap<>();
    for (final SAMReadGroupRecord srgr : header.getReadGroups()) {
      final MachineType mt = ReadGroupUtils.platformToMachineType(srgr, true);
      if (mt == null) {
        throw new NoTalkbackSlimException("Read group with platform specified required");
      } else if (mt.orientation() == null
          || mt != MachineType.ILLUMINA_PE && mt != MachineType.COMPLETE_GENOMICS && mt != MachineType.COMPLETE_GENOMICS_2 //TODO prove other types are ok and this check can be removed
          ) {
        throw new NoTalkbackSlimException("Unsupported platform: " + srgr.getPlatform());
      }
      final String rgId = getReadGroupLabel(srgr.getReadGroupId());
      rgMachineTypes.put(rgId, mt);
    }

    final ReadGroupState[] states = new ReadGroupState[rgStats.size()];
    int s = 0;
    int maxRadius = 0;
    for (final Map.Entry<String, ReadGroupStats> entry : rgStats.entrySet()) {
      final String key = entry.getKey();
      final MachineType machineType = rgMachineTypes.get(key);
      if (machineType == null) {
        throw new NoTalkbackSlimException("Read group " + entry.getKey() + " not contained in SAM header.");
      }
      final ReadGroupState s1 = new ReadGroupState(entry.getValue(), machineType, mCorrections);
      maxRadius = Math.max(Math.max(maxRadius, s1.hi()), -s1.lo());
      states[s++] = s1;
      Diagnostic.userLog("Created read group state for read group: " + s1.stats().id() + " with machine type " + machineType);
      mReadGroupStates.put(entry.getKey(), s1);
    }
    mGlobalRadius = maxRadius + mDefaultStepSize;
    Diagnostic.developerLog("Setting global radius to " + mGlobalRadius);

    if (mParams.outputSimple()) {
      final int lo = -(mParams.binSize() / 2);
      final int hi = (mParams.binSize() + 1) / 2;
      assert hi - lo == mParams.binSize();
      @SuppressWarnings("unchecked")
      final ArrayList<Signal>[] simpleSubSignals = (ArrayList<Signal>[]) new ArrayList<?>[9];
      for (int i = 0; i < simpleSubSignals.length; ++i) {
        simpleSubSignals[i] = new ArrayList<>();
      }
      for (final ReadGroupState rgs : states) {
        simpleSubSignals[0].add(new SignalCount(rgs.properLeftArm(false), lo, hi, "proper-left"));
        simpleSubSignals[1].add(new SignalCount(rgs.discordantLeftArm(false), lo, hi, "discordant-left"));
        simpleSubSignals[2].add(new SignalCount(rgs.unmatedLeftArm(false), lo, hi, "unmated-left"));
        simpleSubSignals[3].add(new SignalCount(rgs.properRightArm(false), lo, hi, "proper-right"));
        simpleSubSignals[4].add(new SignalCount(rgs.discordantRightArm(false), lo, hi, "discordant-right"));
        simpleSubSignals[5].add(new SignalCount(rgs.unmatedRightArm(false), lo, hi, "unmated-right"));
        simpleSubSignals[6].add(new SignalCount(rgs.notPaired(), lo, hi, "not-paired"));
        simpleSubSignals[7].add(new SignalCount(rgs.unique(), lo, hi, "unique"));
        simpleSubSignals[8].add(new SignalCount(rgs.ambiguous(), lo, hi, "ambiguous"));
      }
      mSimpleSignals = new Signal[] {
          new SignalSum("proper-left",      getSignals(simpleSubSignals[0])),
          new SignalSum("discordant-left",  getSignals(simpleSubSignals[1])),
          new SignalSum("unmated-left",     getSignals(simpleSubSignals[2])),
          new SignalSum("proper-right",     getSignals(simpleSubSignals[3])),
          new SignalSum("discordant-right", getSignals(simpleSubSignals[4])),
          new SignalSum("unmated-right",    getSignals(simpleSubSignals[5])),
          new SignalSum("not-paired",       getSignals(simpleSubSignals[6])),
          new SignalSum("unique",           getSignals(simpleSubSignals[7])),
          new SignalSum("ambiguous",        getSignals(simpleSubSignals[8])),
      };
    } else {
      mSimpleSignals = null;
    }


    final ArrayList<Signal> ps = new ArrayList<>();

    // Normal diploid type coverage signals
    ps.add(new NormalBayesianSignal(2).makeSignal(states, false, "normal"));

    // Homozygous models
    ps.add(new NormalBayesianSignal(4).makeSignal(states, false, "duplicate"));
    ps.add(new NormalBayesianSignal(0).makeSignal(states, false, "delete"));
    ps.add(new DeleteBoundaryBayesianSignal().makeSignal(states, false, "delete-left"));
    ps.add(new DeleteBoundaryBayesianSignal().makeSignal(states, true, "delete-right"));
    ps.add(new DuplicateDonorBayesianSignal().makeSignal(states, false, "duplicate-left"));
    ps.add(new DuplicateDonorBayesianSignal().makeSignal(states, true, "duplicate-right"));
    ps.add(new BreakpointBayesianSignal().makeSignal(states, false, "breakpoint"));
    ps.add(new NovelInsertionBayesianSignal().makeSignal(states, true, "novel-insertion"));

    if (mParams.heterozygous()) {
      final NormalBayesianSignal diploidNorm = new NormalBayesianSignal(2);
      // Heterozygous models
      ps.add(new NormalBayesianSignal(3).makeSignal(states, false, "duplicate-hetero"));
      ps.add(new NormalBayesianSignal(1).makeSignal(states, false, "delete-hetero"));
      ps.add(new HeterozygousBayesianSignal(diploidNorm, new DeleteBoundaryBayesianSignal()).makeSignal(states, false, "delete-left-hetero"));
      ps.add(new HeterozygousBayesianSignal(diploidNorm, new DeleteBoundaryBayesianSignal()).makeSignal(states, true, "delete-right-hetero"));
      ps.add(new HeterozygousBayesianSignal(diploidNorm, new DuplicateDonorBayesianSignal()).makeSignal(states, false, "duplicate-left-hetero"));
      ps.add(new HeterozygousBayesianSignal(diploidNorm, new DuplicateDonorBayesianSignal()).makeSignal(states, true, "duplicate-right-hetero"));
      ps.add(new HeterozygousBayesianSignal(diploidNorm, new BreakpointBayesianSignal()).makeSignal(states, false, "breakpoint-hetero"));
      ps.add(new HeterozygousBayesianSignal(diploidNorm, new NovelInsertionBayesianSignal()).makeSignal(states, true, "novel-insertion-hetero"));
    }

    mPosteriorSignals = ps.toArray(new Signal[ps.size()]);

    if (mParams.outputSimple()) {
      writeHeader(mSimpleOut, "simple", mSimpleSignals, "n-count");
    }
    writeHeader(mBayesianOut, "bayesian", mPosteriorSignals, "max-index");
  }

  private Signal[] getSignals(ArrayList<Signal> list) {
    return list.toArray(new Signal[list.size()]);
  }

  private void warnRgMismatch(SAMRecord rec) {
    if (mWarnCount < MAX_WARNINGS) {
      Diagnostic.warning("Skipping record with unrecognized read group: " + rec.getSAMString().trim());
      ++mWarnCount;
      if (mWarnCount == MAX_WARNINGS) {
        Diagnostic.warning("Subsequent warnings of this type will not be shown.");
      }
    }
  }

  @Override
  protected void finalPostFlush() throws IOException {
    processLastTemplate();
  }

  @Override
  public int flush(final int start, final int last) throws IOException {
    assert 0 <= start;
    assert start < last;

    // Ensure mTemplateNs is appropriately populated
    if (mTemplatePos < start - mGlobalRadius) {
      mTemplatePos = start - mGlobalRadius;
      mTemplateNs.flushTo(mTemplatePos);
    }

    int backStepFrom = -1;
    for (int i = start; i < last; ++i) {
      while (mTemplatePos < Math.min(i + mGlobalRadius, mTemplateLength)) {
        if (mTemplate[mTemplatePos] == 0) {
          mTemplateNs.increment(mTemplatePos);
        }
        ++mTemplatePos;
      }
      if (i % mCurrentStepSize == 0) {
        final double nCount = mTemplateNs.count(i, -mGlobalRadius, mGlobalRadius);
        final double[] simple;
        if (mParams.outputSimple()) {
          simple = new double[mSimpleSignals.length + 1]; // Add 1 for n-count
          for (int s = 0; s < mSimpleSignals.length; ++s) {
            final Signal sig = mSimpleSignals[s];
            simple[s] = sig.value(i);
          }
          simple[simple.length - 1] = nCount;
        } else {
          simple = null;
        }
        final double[] norm = normalizedBayes(i);
        final int bestHypothesis = maxIndex(norm);
        if (bestHypothesis != mLastHypothesis) {
          if (mCurrentStepSize > mZoomedStepSize) {
            // We need to step backwards and examine more closely
            backStepFrom = i;
            i = Math.max(i - mDefaultStepSize, mLastWrite); //don't jump before 0
            //System.err.println("backStep from:" + backStepFrom + " to: " + i);
            mCurrentStepSize = mZoomedStepSize;
            mLastChange = i;
            continue;
          }
          // We're already zoomed so proceed remembering where we last changed
          mLastChange = i;
        } else if (mLastChange < i - mDefaultStepSize - 1 && mCurrentStepSize < mDefaultStepSize) {
          // we've consumed a big step size without changing the best
          // hypothesis so go back to large jumps
          mCurrentStepSize = mDefaultStepSize;
        }

        final boolean inN = nCount > mGlobalRadius; // 50% N's
        mLastHypothesis = bestHypothesis;
        if (!inN) {
          writeValues(mBayesianOut, i, norm, true);
        }
        mLastWrite = i;
        if (simple != null) {
          writeValues(mSimpleOut, i, simple, false);
        }
      }
      if (i >= backStepFrom) {
        // We don't want to flush things that are ahead of us on a backstep
        flushStates(i - mGlobalRadius);
      }
    }
    return last;
  }

  private static int maxIndex(double[] arr) {
    double max = Double.NEGATIVE_INFINITY;
    int maxIndex = 0;
    for (int i = 0; i < arr.length; ++i) {
      if (arr[i] > max) {
        max = arr[i];
        maxIndex = i;
      }
    }
    return maxIndex;
  }


  private void resetStates(int start) {
    mTemplateNs.reset(mTemplateLength, (mGlobalRadius + 1) * 2);
    for (final ReadGroupState state : mReadGroupStates.values()) {
      state.reset(mTemplateLength, (mGlobalRadius + 1) * 2);
    }
    if (start - mGlobalRadius > 0) {
      flushStates(start - mGlobalRadius);
    }
  }
  private void flushStates(int offset) {
    if (offset >= 0) {
      mTemplateNs.flushTo(offset);
      for (final ReadGroupState state : mReadGroupStates.values()) {
        state.flushTo(offset);
      }
    }
  }


  private void writeValues(OutputStream out, int zeroBasedPosition, double[] values, boolean addMax) throws IOException {
    out.write(mTemplateName.getBytes());
    out.write(TAB.getBytes());
    out.write(("" + (zeroBasedPosition + 1)).getBytes());
    for (final double value : values) {
      out.write(TAB.getBytes());
      out.write(Utils.realFormat(value, 4).getBytes());
    }
    if (addMax) {
      out.write(TAB.getBytes());
      final int maxIndex = maxIndex(values);
      out.write(Integer.toString(maxIndex).getBytes());
      mInterestingRegionExtractor.processValue(zeroBasedPosition + 1, values, maxIndex);
    }
    out.write(StringUtils.LS.getBytes());
  }

  private double[] normalizedBayes(final int position) {
    final double[] posterior = new double[mPosteriorSignals.length];
    for (int i = 0; i < mPosteriorSignals.length; ++i) {
      posterior[i] = mPosteriorSignals[i].value(position);
    }
    if (posterior.length == 1) {
      return posterior; /// Just here for debugninn
    }
    return normalize(posterior);
  }

  static double[] normalize(final double[] logValues) {
    final int length = logValues.length;
    double min = Double.POSITIVE_INFINITY;

    int mini = -1;
    for (int i = 0; i < length; ++i) {
      final double lv = logValues[i];
      assert Double.isFinite(lv);
      if (lv < min) {
        min = lv;
        mini = i;
      }
    }
    assert mini != -1;

    final double[] norm = new double[length];
    double othersum = Double.NEGATIVE_INFINITY;
    double allsum = Double.NEGATIVE_INFINITY;
    for (int i = 0; i < length; ++i) {
      if (i != mini) {
        othersum = VariantUtils.logSumApproximation(othersum, -logValues[i]);
      }
      allsum = VariantUtils.logSumApproximation(allsum, -logValues[i]);
    }

    final double log10 = Math.log(10);
    for (int i = 0; i < length; ++i) {
      final double lv = -logValues[i];
      assert Double.isFinite(lv);
      if (i == mini) {
        norm[i] = (lv - othersum) / log10;
      } else {
        norm[i] = (lv - VariantUtils.logSubtract(allsum, lv)) / log10;
      }
      assert Double.isFinite(norm[i]);
    }
    return norm;
  }

  @Override
  protected boolean processRecord(final SAMRecord rec) throws IOException {
    final int start = rec.getAlignmentStart() - 1; // zero-based, inclusive
    final int end = rec.getAlignmentEnd();  // zero-based, exclusive
    final String refName = rec.getReferenceName();
    if (!refName.equals(mTemplateName)) {
      if (mTemplateName != null) {
        flush(mPreviousStart, mTemplateLength);
        //processLastTemplate();
      }
      final Long templateId = mTemplateNameMap.get(refName);
      if (templateId == null) {
        throw new NoTalkbackSlimException("Wrong reference sequences supplied for the given mappings.");
      }
      final int templateLength = getSequenceLength(mGenomeSequences, mTemplateNameMap, refName);
      final String restrictionTemplate = mParams.filterParams().restrictionTemplate();
      if ((restrictionTemplate != null) && restrictionTemplate.equals(refName) && mParams.filterParams().restrictionStart() != -1) {
        mTemplateLength = mParams.filterParams().restrictionEnd();
        mTemplateMin = mParams.filterParams().restrictionStart();
      } else {
        mTemplateLength = templateLength;
        mTemplateMin = 0;
      }
      mInterestingRegionExtractor.setTemplate(refName, templateLength);
      mTemplate = mGenomeSequences.read(templateId);
      mTemplateName = refName;
      mTemplatePos = mTemplateMin;
      mLastWrite = mTemplateMin;
      mPreviousStart = mTemplateMin;
      resetStates(mTemplateMin);
    } else {
      if (start < mPreviousSam) {
        throw new NoTalkbackSlimException(ErrorType.SAM_NOT_SORTED);
      }
    }
    mPreviousSam = start;
    if (end < mTemplateMin || start > mTemplateLength) {
      // does not intersect template at all.
      return false;
    } else if (start >= end) {
      // Invalid range
      return false;
    }

    if (start - mGlobalRadius > mPreviousStart) {
      flush(mPreviousStart, start - mGlobalRadius);
      mPreviousStart = start - mGlobalRadius;
    }

    final ReadGroupState rgstate = getState(rec);
    if (rgstate == null) {
      warnRgMismatch(rec);
      return false;
    } else {
      return rgstate.update(rec);
    }
  }

  private ReadGroupState getState(SAMRecord rec) {
    return mReadGroupStates.get(getReadGroupLabel(rec));
  }

  private void processLastTemplate() throws IOException {
    mInterestingRegionExtractor.setTemplate(null, -1);
  }

  private void writeHeader(OutputStream out, String outputType, Signal[] columns, String extra) throws IOException {
    out.write(("#Version " + Environment.getVersion() + ", SV " + outputType + " output " + SV_OUTPUT_VERSION + StringUtils.LS).getBytes());
    if (CommandLine.getCommandLine() != null) {
      out.write(("#CL" + TAB + CommandLine.getCommandLine() + StringUtils.LS).getBytes());
    }
    out.write(("#RUN-ID" + TAB + CommandLine.getRunId() + StringUtils.LS).getBytes());
    out.write(("#template-name" + TAB + "position").getBytes());
    for (final Signal sig : columns) {
      out.write((TAB + sig.columnLabel()).getBytes());
    }
    if (extra != null) {
      out.write((TAB + extra).getBytes());
    }
    out.write(StringUtils.LS.getBytes());
  }

  @Override
  protected void exec() throws IOException {
    mSimpleOut = mParams.outputSimple() ? mParams.outStream(SvToolParams.NAME_SIMPLE) : null;
    try {
      mBayesianOut = mParams.bayesianStream();
      try {
        try {
          mInterestingRegionExtractor = new SvInterestingRegionExtractor(mParams);
          try {
            super.exec();
          } finally {
            mInterestingRegionExtractor.close();
          }
        } finally {
          mGenomeSequences.close();
        }
      } finally {
        mBayesianOut.close();
      }
    } finally {
      if (mParams.outputSimple()) {
        mSimpleOut.close();
      }
    }
  }
}
