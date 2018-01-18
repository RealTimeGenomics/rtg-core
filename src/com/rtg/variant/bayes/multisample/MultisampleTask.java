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
package com.rtg.variant.bayes.multisample;

import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import com.rtg.bed.BedUtils;
import com.rtg.launcher.ParamsTask;
import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.reference.Ploidy;
import com.rtg.reference.SexMemo;
import com.rtg.sam.CircularBufferMultifileSinglePassReaderWindowSync;
import com.rtg.sam.ReaderRecord;
import com.rtg.sam.ReaderWindow;
import com.rtg.sam.SamReadingContext;
import com.rtg.sam.SamUtils;
import com.rtg.sam.ThreadedMultifileIteratorWrapper;
import com.rtg.scheduler.EventList;
import com.rtg.scheduler.Executor;
import com.rtg.scheduler.ExecutorRandom;
import com.rtg.scheduler.ExecutorSequential;
import com.rtg.scheduler.ExecutorThreaded;
import com.rtg.scheduler.Job;
import com.rtg.scheduler.JobFactory;
import com.rtg.scheduler.JobStatistics;
import com.rtg.scheduler.Result;
import com.rtg.scheduler.Scheduler;
import com.rtg.scheduler.SchedulerSynchronized;
import com.rtg.usage.UsageMetric;
import com.rtg.util.SingletonPopulatorFactory;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.ParallelProgress;
import com.rtg.util.diagnostic.WarningType;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.util.intervals.RangeList;
import com.rtg.util.intervals.ReferenceRegions;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.util.intervals.StatusInterval;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.Variant;
import com.rtg.variant.Variant.VariantFilter;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantAlignmentRecordPopulator;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.avr.ModelFactory;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.complex.AligningDecomposer;
import com.rtg.variant.bayes.complex.Decomposer;
import com.rtg.variant.bayes.complex.DoNothingDecomposer;
import com.rtg.variant.bayes.complex.IonTorrentCallFilter;
import com.rtg.variant.bayes.complex.SimpleDecomposer;
import com.rtg.variant.bayes.complex.Trimmer;
import com.rtg.variant.bayes.multisample.multithread.DependenciesMultiSample;
import com.rtg.variant.bayes.multisample.multithread.EventListMultiSample;
import com.rtg.variant.bayes.multisample.multithread.JobIdMultisample;
import com.rtg.variant.bayes.multisample.multithread.MultisampleStatistics;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.format.VariantOutputVcfFormatter;
import com.rtg.variant.util.VariantUtils;
import com.rtg.vcf.StatisticsVcfWriter;
import com.rtg.vcf.VariantStatistics;
import com.rtg.vcf.VcfAnnotator;
import com.rtg.vcf.VcfFilter;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;
import com.rtg.vcf.VcfWriter;
import com.rtg.vcf.VcfWriterFactory;
import com.rtg.vcf.header.VcfHeader;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

/**
 * @param <V> type of variant statistics
 */
public class MultisampleTask<V extends VariantStatistics> extends ParamsTask<VariantParams, V> {

  private static final int ION_TORRENT_HYPER_COMPLEX = 21;

  private static final int MIN_CALLS_FOR_COVERAGE_WARNING = 50; // Only warn for non-trivial datasets
  private static final double COVERAGE_WARNING_THRESHOLD = 1.0;

  private final JobStatistics<JobIdMultisample> mJobStatistics = new MultisampleStatistics();
  private final OutputStream mBedOut;
  private final SexMemo mSexMemo;
  private final List<VcfAnnotator> mAnnotators = new ArrayList<>();
  private final List<VcfFilter> mFilters = new ArrayList<>();
  protected final SequencesReader mReferenceSequences;
  private AbstractJointCallerConfiguration mConfig;
  private VcfWriter mOut;
  private VariantOutputVcfFormatter mFormatter;
  private VcfHeader mVcfHeader;
  private ThreadedMultifileIteratorWrapper<VariantAlignmentRecord> mWrapper;
  private List<SAMSequenceRecord> mSequences;
  private ReferenceRegions mBedFilterRegions;
  private ParallelProgress mPP = null;
  private Decomposer mDecomposer = null;

  private final JointCallerConfigurator<V> mConfigurator;

  /**
   * @param params command line parameters.
   * @param configurator the joint caller configuration creator
   * @param defaultOutput where the calling results are written.
   * @param statistics the variant statistics collector
   * @param usageMetric count of the number of records read.
   * @throws IOException when writing output.
   */
  public MultisampleTask(VariantParams params, JointCallerConfigurator<V> configurator, OutputStream defaultOutput, V statistics, UsageMetric usageMetric) throws IOException {
    super(params, defaultOutput, statistics, usageMetric);
    if (params.genome() == null) {
      throw new NoTalkbackSlimException(ErrorType.READING_ERROR, "Problem reading");
    }
    mReferenceSequences = params.genome().reader();
    mBedOut = params.bedStream();
    mSexMemo = Utils.createSexMemo(mParams);
    mConfigurator = configurator;

    statistics.setReference(mReferenceSequences.getReadMe() == null ? mReferenceSequences.path().getPath() : mReferenceSequences.getReadMe());
    Diagnostic.developerLog("Genome priors:" + StringUtils.LS + VariantUtils.toGenomePriorProperties(params.genomePriors()));
  }

  private Variant stepIndels(String refName, IndividualSampleProcessor<?>[] ssProcessors, int pos) {
    Variant v = null;
    for (IndividualSampleProcessor<?> processor : ssProcessors) {
      final Variant current = processor.indelOutput(refName, pos, mParams);
      if (v == null || (current != null && current.getLocus().getLength() > v.getLocus().getLength())) {
        v = current;
      }
    }
    return v;
  }

  private static final byte OVERFLOW = 1;
  private static final byte SKIP = 2;

  private int processNtPositions(List<Variant> calls, MultisampleJointCaller jointCaller, ChunkInfo chunkInfo, byte[] template, ReaderWindow<VariantAlignmentRecord> tribble, int start, int end) throws IOException {

    int maxReadLen = 0;
    List<RangeList.RangeData<String>> ranges = null;
    int rangeIndex = 0;
    boolean skipWholeChunk = false;
    if (mWrapper.context().hasRegions()) {
      final RangeList<String> currentRangeList = mWrapper.getCurrentRangeList();
      ranges = currentRangeList.getFullRangeList();
      rangeIndex = currentRangeList.findFullRangeIndex(start);
      assert rangeIndex < ranges.size();
      final RangeList.RangeData<String> range = ranges.get(rangeIndex);
      skipWholeChunk = range.getMeta() == null && range.isInRange(end);
    }
    if (skipWholeChunk) {
      tribble.advanceBuffer(end);
      //Diagnostic.developerLog("Entirely skipping chunk " + chunkInfo.toString() + " due to region " + range.toString());
    } else { // Don't bother pulling records if our entire chunk is contained in a don't-call block

      final Iterator<VariantAlignmentRecord> it = tribble.recordsOverlap(start, end);
      final String refName = chunkInfo.refName();
      final IndividualSampleProcessor<?>[] ssProcessors = mConfig.getIndividualSampleProcessors(refName, template, start, end);

      final StatusInterval statusInterval = new StatusInterval(start, end);
      int overflowEnd = 0; // keeps track of current overflow status

      boolean sawValidRecord = false;
      while (it.hasNext()) {
        final VariantAlignmentRecord rec = it.next();
        //System.err.println(rec);
        final int recStart = rec.getStart();
        final int recEnd = recStart + rec.getLength();
        if (rec.isOverflow()) {
          // Update overflow end position to be largest seen so far (allows for overlapping overflow records)
          statusInterval.add(recStart, recEnd, OVERFLOW);
          overflowEnd = Math.max(overflowEnd, recEnd);
          //Diagnostic.developerLog("Suppressing calls in [" + recStart + "," + overflowEnd + ")");
        } else if (recEnd > overflowEnd) {
          maxReadLen = Math.max(maxReadLen, rec.getLength());
          final IndividualSampleProcessor<?> individualSampleProcessor = ssProcessors[mConfig.numberOfGenomes() == 1 ? 0 : rec.getGenome()];
          if (!individualSampleProcessor.processAlignmentRecord(rec)) {
            addInvalidRecord(rec);
          } else {
            sawValidRecord = true;
          }
        }
      }
      if (sawValidRecord) { // Only bother to step and call in models if we saw some evidence

        // Add in ranges to the status interval as positions to skip
        if (ranges != null) {
          addRangeStatuses(statusInterval, ranges, rangeIndex, end);
        }

        final List<ModelInterface<?>> models = new ArrayList<>(ssProcessors.length);
        for (int pos = start; pos < end; ) {
          if (statusInterval.contains(pos)) {
            final byte status = statusInterval.get(pos);
            final int oldpos = pos;
            do {
              for (final IndividualSampleProcessor<?> ssProcessor1 : ssProcessors) {
                ssProcessor1.step(pos);
              }
              for (final IndividualSampleProcessor<?> ssProcessor : ssProcessors) {
                ssProcessor.stepIndel(pos);
              }
            } while (++pos < end && statusInterval.get(pos) == status);
            if (status == OVERFLOW) {
              final VariantLocus overflowLocus = new VariantLocus(refName, oldpos, pos);
              final Variant v = new Variant(overflowLocus);
              v.setOverflow();
              calls.add(v);
            } // Otherwise status == SKIP, do nothing
          } else {
            final Variant v = stepIndels(refName, ssProcessors, pos);
            if (v != null) {
              calls.add(v);
            }
            models.clear();
            for (final IndividualSampleProcessor<?> ssProcessor : ssProcessors) {
              models.add(ssProcessor.step(pos));
            }

            int refHyp = Integer.MIN_VALUE;
            for (ModelInterface<?> m : models) {
              if (m.hypotheses().size() > 0) {
                refHyp = m.reference(); // todo This is still a bit icky. We are assuming that if you have hypotheses all references are the same (maybe not true in cancer)
                break; // todo this break essential for allele based cancer calling -- check it doesn't break anything else
              }
            }
            if (refHyp != Integer.MIN_VALUE) {
              final HaploidDiploidHypotheses<HypothesesPrior<Description>> hypotheses = mConfig.getSnpHypotheses(refHyp, refName, pos);
              final Variant variant;
              final boolean onlyRefCoverage = Utils.hasOnlyRefCoverage(models);
              if (mParams.callLevel() != VariantOutputLevel.ALL && onlyRefCoverage) {
                variant = null;
              } else {
                for (ModelInterface<?> model : models) {
                  model.freeze();
                }
                variant = jointCaller.makeCall(refName, pos, pos + 1, template, models, hypotheses);
              }

              if (variant != null) {
                calls.add(variant);
              }
            }
            ++pos;
          }
        }
      }
    }
    mPP.updateProgress(chunkInfo.percent(end));
    return maxReadLen;
  }

  // Sets the status of any positions within the interval to SKIP if they are contained within a no-call range entry (one without metadata)
  private static void addRangeStatuses(StatusInterval statusInterval, List<RangeList.RangeData<String>> ranges, int startIndex, int endPos) {
    for (int rangeIndex = startIndex; rangeIndex < ranges.size(); ++rangeIndex) {
      final RangeList.RangeData<String> range = ranges.get(rangeIndex);
      if (range.getStart() > endPos) {
        break;
      }
      if (range.getMeta() == null) {
        statusInterval.add(range.getStart(), range.getEnd(), SKIP);
      }
    }
  }

  private void addInvalidRecord(ReaderRecord<?> record) {
    ++mInvalidRecords;
    if (mInvalidRecords <= 5) {
      Diagnostic.warning(WarningType.SAM_BAD_FORMAT_WARNING1, record.toString());
      Diagnostic.userLog("Invalid record: " + record);
    }
  }

  private void updateCounts(final ComplexCaller caller) {
    synchronized (mExcessiveCoverageLock) {
      mExcessiveCoverageCount += caller.getExcessiveCoverageCount();
    }
    synchronized (mExcessiveHypothesesLock) {
      mExcessiveHypothesesCount += caller.getExcessiveHypothesesCount();
    }
    synchronized (mNoHypothesesLock) {
      mNoHypothesesCount += caller.getNoHypothesesCount();
    }
  }

  private long mInvalidRecords = 0;
  private long mTossedRecords = 0;
  private long mExcessiveCoverageCount = 0;
  private long mExcessiveHypothesesCount = 0;
  private long mNoHypothesesCount = 0;

  private final Object mExcessiveCoverageLock = new Object();
  private final Object mExcessiveHypothesesLock = new Object();
  private final Object mNoHypothesesLock = new Object();

  private class JobFactoryMultiSample extends IntegralAbstract implements JobFactory<JobIdMultisample>, AutoCloseable {

    private final ChunkInfo mInfo;
    private final CircularBufferMultifileSinglePassReaderWindowSync<VariantAlignmentRecord> mBuffer;
    private final String mRefName;
    private final byte[] mRefNts;
    private final MultisampleJointCaller mJointCaller;
    private final BedComplexitiesWriter mBed;

    /** Minimum position on reference, either 0 or supplied by restriction */
    private final int mMinimumPosition;

    JobFactoryMultiSample(final ChunkInfo info, final String refName, final byte[] refNts) throws IOException {
      mInfo = info;
      String[] genomeNames = mConfig.getGenomeNames();
      if (genomeNames.length == 1) {
        genomeNames = new String[] {}; // Special case for singleton caller, map all records to 0
      }
      final int depth = mParams.maxCoverageBypass().thresholdTotal(refName);
      final VariantAlignmentRecordPopulator pop = new VariantAlignmentRecordPopulator(MultisampleUtils.chooser(mParams), mParams.minBaseQuality(), genomeNames);
      final RegionRestriction restriction = new RegionRestriction(refName, info.start(), info.end());

      if (restriction.getStart() < 0) {
        mMinimumPosition = 0;
      } else {
        mMinimumPosition = restriction.getStart();
      }
      mBuffer = new CircularBufferMultifileSinglePassReaderWindowSync<>(mWrapper, pop, mParams.uberHeader().getSequenceIndex(refName), restriction.getStart(), depth);
      mRefName = refName;
      mRefNts = refNts;
      mJointCaller = mConfig.getJointCaller();
      mBed = new BedComplexitiesWriter(mBedOut, refName, info.start());
    }

    @Override
    public void close() throws IOException {
      mTossedRecords += mBuffer.getTossedRecordCount();
      mBed.finish();
      mBuffer.close();
      mPP.close();
      mJointCaller.endOfSequence();
    }

    @Override
    public boolean integrity() {
      Exam.assertTrue(mInfo.numberChunks() > 0);
      return true;
    }

    @Override
    public Job<JobIdMultisample> job(final JobIdMultisample id, final Result[] arguments) {
      switch (id.type()) {
        case INCR:
          return new IncrJob(id, arguments);
        case DANGLING:
          return new DanglingJob(id, arguments);
        case BED:
          return new BedJob(id, arguments);
        case COMPLEX:
          return new ComplexJob(id, arguments);
        case FLUSH:
          return new FlushJob(id, arguments);
        case FILTER:
          return new FilterJob(id, arguments);
        case OUT:
          return new OutJob(id, arguments);
        default:
          throw new RuntimeException();
      }
    }

    private abstract class MultisampleJob extends Job<JobIdMultisample> {
      protected final Result[] mArguments;
      private final int mTimeOffset;
      MultisampleJob(JobIdMultisample id, Result[] arguments, int timeOffset) {
        super(id);
        mArguments = arguments;
        mTimeOffset = timeOffset;
      }

      @Override
      public String toString() {
        // The positions reported here are only a proxy for the true positions, since for example, the DANGLING
        // job can make adjustments to the positions processed.  The idea is to try and make the position
        // reported here correspond to the original INCR job used to generate the inputs of this chunk.
        final int delta = Math.max((id().time() - mTimeOffset) * mInfo.chunkSize(), 0);
        final int start = Math.min(delta, mInfo.end());
        final int end = Math.min(start + mInfo.chunkSize(), mInfo.end());
        return super.toString() + " " + mRefName + ":" + start + "-" + end;
      }

      @SuppressWarnings("unchecked")
      protected List<Variant> getList(Object o) {
        return (List<Variant>) o;
      }

    }

    private class IncrJob extends MultisampleJob {
      IncrJob(JobIdMultisample id, Result[] arguments) {
        super(id, arguments, 0);
      }

      @Override
      public Result run() throws IOException {
        final int start = id().time() * mInfo.chunkSize() + mInfo.start();
        final int end = Math.min(start + mInfo.chunkSize(), mInfo.end());
        List<Variant> calls = new ArrayList<>();
        final int maxReadLen = processNtPositions(calls, mJointCaller, mInfo, mRefNts, mBuffer, start, end);
        final boolean simpleRepeats = mParams.simpleRepeatExtension() && !mParams.ionTorrent();
        final RegionRestriction forcedComplexRegion = mParams.forceComplexRegion();
        if (forcedComplexRegion != null) {
          if (forcedComplexRegion.getSequenceName().equals(mRefName) && start <= forcedComplexRegion.getStart() && end > forcedComplexRegion.getStart()) {
            final Variant v = new Variant(new VariantLocus(mRefName, forcedComplexRegion.getStart(), forcedComplexRegion.getEnd()));
            v.setInteresting();
            v.setIndel(1); //forces complex
            calls = OutputUtils.merge(calls, Collections.singletonList(v));
          }
        }
        final Complexities cx = new Complexities(calls, mRefName, start, end, mParams.interestingSeparation(), mParams.hyperComplexLength(), mRefNts, simpleRepeats, mConfig.getSiteSpecificPriors());
        return new Result(cx, maxReadLen);
      }
    }

    private class DanglingJob extends MultisampleJob {
      DanglingJob(JobIdMultisample id, Result[] arguments) {
        super(id, arguments, 1);
      }

      @Override
      public Result run() {
        final Complexities last = mArguments[0] == null ? null : (Complexities) mArguments[0].result(0);
        final Complexities cx = mArguments[1] == null ? null : (Complexities) mArguments[1].result(0);
        Complexities.fixDangling(last, cx);
        //passing through complexities from previous incr time step to current complex time step
        return new Result(last);
      }
    }

    private class BedJob extends MultisampleJob {
      BedJob(JobIdMultisample id, Result[] arguments) {
        super(id, arguments, 1);
      }

      @Override
      public Result run() throws IOException {
        final Complexities last = (Complexities) mArguments[1].result(0);
        if (last != null) {
          mBed.write(last);
        }
        return new Result();
      }
    }

    private class ComplexJob extends MultisampleJob {
      ComplexJob(JobIdMultisample id, Result[] arguments) {
        super(id, arguments, 1);
      }

      @Override
      public Result run() throws IOException {
        final Complexities complexRegions = (Complexities) mArguments[0].result(0);
        if (complexRegions != null) {
          assert complexRegions.isFixed();
          final List<Variant> calls;
          if (mParams.noComplexCalls()) {
            calls = new ArrayList<>();
            for (final Variant v : complexRegions.getOriginalCalls()) {
              if (!v.isIndel() && !v.isSoftClip() && !v.isOverflow() && v.getLocus().getStart() >= complexRegions.startOfChunk() && v.getLocus().getStart() < complexRegions.endOfChunk()) {
                calls.add(v);
              }
            }
          } else {
            final ComplexCaller caller = new ComplexCaller(mParams, mConfig /*, mWrapper.getCurrentRangeList() */);
            final List<Variant> complexCalls = caller.makeComplexCalls(complexRegions, mBuffer, mRefNts, mRefName);
            final List<Variant> nonComplexCalls = OutputUtils.nonShadowed(complexRegions.getOriginalCalls(), complexRegions);
            calls = OutputUtils.merge(nonComplexCalls, complexCalls);
            updateCounts(caller);
          }
          return new Result(calls, complexRegions);
        } else {
          return new Result(null, null);
        }
      }
    }

    private class FlushJob extends MultisampleJob {
      FlushJob(JobIdMultisample id, Result[] arguments) {
        super(id, arguments, 1);
      }

      @Override
      protected Result run() throws IOException {
        final Complexities last = mArguments[0] == null ? null : (Complexities) mArguments[0].result(0);
        //operating 1 time step behind increment job. So operating on 'last'
        if (last != null) {
          if (mParams.expandComplexReadQueries()) {
            final Complexities current = mArguments[1] == null ? null : (Complexities) mArguments[1].result(1);
            //complex calls may have are will be made +/- 1 outside of chunk boundaries
            final int flushStart = Math.max(mMinimumPosition, last.startOfChunk() - 1);
            //don't flush last base so it can be used by next complex region
            final int flushEnd = current != null ? last.endOfChunk() - 1 : last.endOfChunk();
            mBuffer.flush(flushStart, flushEnd);
          } else {
            mBuffer.flush(last.startOfChunk(), last.endOfChunk());
          }
        }
        return new Result();
      }
    }

    private class FilterJob extends MultisampleJob {
      FilterJob(JobIdMultisample id, Result[] arguments) {
        super(id, arguments, 1);
      }

      @Override
      public Result run() {
        final List<Variant> prevLastCall = mArguments[2] == null ? null : getList(mArguments[2].result(1));
        final List<Variant> equivFiltered;
        final List<Variant> lastCalls;
        if (mArguments[0] == null && mArguments[1] == null) { // Final filter job to flush out previous last call (if present)
          equivFiltered = prevLastCall;
          lastCalls = null;
        } else {
          final Integer maxReadLen = (Integer) mArguments[0].result(1);
          final List<Variant> initialCalls = getList(mArguments[1].result(0));
          final List<Variant> split = trimSplit(initialCalls);
          final List<Variant> filtered = locusAndIonTorrentFilters(split, mRefNts);
          final EquivalentFilter filter = new EquivalentFilter(mRefNts, prevLastCall);
          equivFiltered = filter.filter(filtered, maxReadLen);
          lastCalls = filter.lastCall(); // Remember last calls for checking equivalence across chunks
        }
        if (equivFiltered != null) {
          final List<VcfRecord> vcfRecords = new ArrayList<>(equivFiltered.size());
          for (final Variant v : equivFiltered) {
            final VcfRecord record = mFormatter.makeVcfRecord(v);
            for (final VcfAnnotator annot : mAnnotators) {
              annot.annotate(record);
            }
            vcfRecords.add(record);
          }
          return new Result(vcfRecords, lastCalls);
        } else {
          return new Result(null, lastCalls);
        }
      }
    }

    private class OutJob extends MultisampleJob {
      OutJob(JobIdMultisample id, Result[] arguments) {
        super(id, arguments, 2);
      }

      @Override
      public Result run() throws IOException {
        final List<VcfRecord> filtered = getVcfList(mArguments[1].result(0));
        if (filtered != null) {
          for (final VcfRecord record : filtered) {
            boolean keep = true;
            for (final VcfFilter filter : mFilters) {
              if (!filter.accept(record)) {
                keep = false;
                break;
              }
            }
            if (keep) {
              mOut.write(record);
            }
          }
        }
        if (mBuffer.finishedTo() < Math.min(id().time() * mInfo.chunkSize() + mInfo.start() - 1, mInfo.end())) { //flushing should be keeping up with output
          throw new RuntimeException("Failed to flush chunk: " + mBuffer.finishedTo() + " : " + id().time() * mInfo.chunkSize());
        }
        return null;
      }
    }
  }

  private List<Variant> trimSplit(List<Variant> merged) {
    final List<Variant> calls = new ArrayList<>();
    for (Variant variant : merged) {
      if (variant.isOrdinaryCall() && variant.isComplexScored()) {
        calls.addAll(mDecomposer.decompose(variant));
      } else {
        calls.add(variant);
      }
    }
    return calls;
  }


  List<Variant> locusAndIonTorrentFilters(List<Variant> variants, byte[] template) {
    if (variants == null) {
      return null;
    }
    final List<Variant> result = new ArrayList<>();
    for (Variant v : variants) {
      if (!mParams.outputNonSnps() && !v.isSnp()) {
        // User doesn't want non-snps
        continue;
      }
      if (v.getLocus().getStart() == 0 && VariantOutputVcfFormatter.includePreviousNt(v)) {
        // At the start of a reference sequence we currently append an N instead of anchoring on the following base -- for now don't output these records
        continue;
      }
      result.add(v);
      if (IonTorrentCallFilter.ionTorrentFilter(template, v, mParams.ionTorrent())) {
        v.addFilter(VariantFilter.IONTORRENT);
      }
      if (mBedFilterRegions != null && !mBedFilterRegions.overlapped(v.getLocus())) {
        v.addFilter(VariantFilter.BED_REGION);
      }
    }
    return result;
  }

  private void processAnEntireSequence(final String refName, final byte[] refNts) throws IOException {
    mPP = new ParallelProgress(refName);
    final Ploidy ploidy = mSexMemo.getRealPloidy(mParams.sex(), refName);
    if (!mConfig.handlesPloidy(ploidy)) {
      Diagnostic.userLog(ploidy + " sequence " + refName + " not supported in this caller");
      return;
    }
    Diagnostic.userLog("Sequence " + refName + " filter on maximum per-sample coverage is " + mParams.maxCoverageFilter().thresholdSingle(refName));
    Diagnostic.userLog("Sequence " + refName + " extreme coverage bypass level is "
                       + mParams.maxCoverageBypass().thresholdTotal(refName));

    final List<RangeList.RangeData<String>> ranges = mWrapper.getCurrentRangeList().getRangeList();
    assert !ranges.isEmpty();
    final int startPos = ranges.get(0).getStart();
    if (startPos >= refNts.length) {
      throw new NoTalkbackSlimException("Desired start position for sequence " + refName + " (" + startPos + ") is greater than available reference SDF sequence length (" + refNts.length + ")");
    }
    int endPos = ranges.get(ranges.size() - 1).getEnd();
    if (endPos > refNts.length) {
      Diagnostic.warning("Sequence length disparity between reference SDF and SAM headers for sequence " + refName + ". Clipping end position to available SDF sequence length (" + refNts.length + ")");
      endPos = refNts.length;
    }
    final ChunkInfo info = new ChunkInfo(refNts.length, refName, mParams.chunkSize(), startPos, endPos, mParams.execThreads(), mParams.maxReadLength());
    final DependenciesMultiSample depen = new DependenciesMultiSample(info.numberChunks());
    try (final JobFactoryMultiSample jobFac = new JobFactoryMultiSample(info, refName, refNts)) {
      final EventList<JobIdMultisample> eventList = new EventListMultiSample<>();
      final SchedulerSynchronized<JobIdMultisample> sched = new SchedulerSynchronized<>(depen, jobFac, eventList, null, mJobStatistics, mParams.threadingLookAhead());
      //final Scheduler<JobIdMultisample> sched = new SchedulerSynchronized<>(depen, jobFac, eventList, System.err, mJobStatistics, mParams.threadingLookAhead());
      final String msg = "Processing " + refName;
      final Executor<JobIdMultisample> exec = createExecutor(sched, msg, mParams);
      exec.run();
      sched.dumpStarvation();
      assert eventList.next(sched.lookAhead()) == null;
      assert sched.lookAhead().total() == 0;
    }
  }

  private static Executor<JobIdMultisample> createExecutor(final Scheduler<JobIdMultisample> sched, final String msg, VariantParams params) {
    final Executor<JobIdMultisample> exec;
    final String paraMsg;
    switch (params.threadingEnvironment()) {
      case SINGLE:
        exec = new ExecutorSequential<>(sched);
        paraMsg = "single";
        break;
      case RANDOM:
        final Long seed = params.threadingEnvironmentSeed();
        exec = new ExecutorRandom<>(sched, params.execThreads(), seed);
        paraMsg = "random" + (seed == null ? "" : " " + seed);
        break;
      case PARALLEL:
        //exec = new ExecutorSequential(sched); //handy for testing
        //exec = new ExecutorRandom(sched, 10, 42L);
        exec = new ExecutorThreaded<>(sched, params.execThreads());
        paraMsg = "parallel " + params.execThreads();
        break;
      default:
        throw new RuntimeException();
    }
    Diagnostic.developerLog(msg + " " + paraMsg);
    return exec;
  }

  private void logRecordCounts() {
    mUsageMetric.incrementMetric(mWrapper.getTotalNucleotides());
    final long invalidRecords = mInvalidRecords + mWrapper.getInvalidRecordsCount();

    Diagnostic.userLog(mWrapper.getTotalRecordsCount() + " total input alignments");
    final long filteredRecords = mWrapper.getFilteredRecordsCount();
    Diagnostic.userLog(filteredRecords + " alignments skipped due to input filtering criteria");
    final String invalidRecordsWarning = invalidRecords + " alignments skipped because of SAM format problems.";
    if (invalidRecords > 0) {
      Diagnostic.warning(invalidRecordsWarning);
    } else {
      Diagnostic.userLog(invalidRecordsWarning);
    }
    final long numDeduped = mWrapper.getDuplicateRecordsCount();
    Diagnostic.userLog(numDeduped + " alignments skipped due to duplicate detection");
    Diagnostic.userLog(mTossedRecords + " alignments skipped in extreme coverage regions");
    final long processed = mWrapper.getTotalRecordsCount() - invalidRecords - filteredRecords - numDeduped - mTossedRecords;
    Diagnostic.userLog(processed + " alignments processed");
  }



  void init() throws IOException {
    boolean detectedIonTorrent = false;
    for (final SAMReadGroupRecord readGroup : mParams.uberHeader().getReadGroups()) {
      if (MachineType.IONTORRENT.compatiblePlatform(readGroup.getPlatform())) {
        detectedIonTorrent = true; //see bug 1447
      }
    }
    if (detectedIonTorrent) {
      Diagnostic.developerLog("Detected IONTORRENT platform");
      final VariantParamsBuilder vpb = mParams.cloneBuilder();
      vpb.ionTorrent(true);
      if (mParams.hyperComplexLength() == Integer.MAX_VALUE) {
        vpb.hyperComplexLength(ION_TORRENT_HYPER_COMPLEX);
      }
      vpb.pruneHypotheses(true);
      mParams = vpb.create();
      Diagnostic.developerLog("New params:" + mParams);
    }

    mConfig = mConfigurator.getConfig(mParams, mStatistics);
    final VariantAlleleTrigger variantAlleleTrigger = new VariantAlleleTrigger(mParams.minVariantAllelicDepth(), mParams.minVariantAllelicFraction());
    final DecomposerType trimSplitType = mParams.trimSplit();
    if (trimSplitType == DecomposerType.NONE || mParams.callLevel() == VariantOutputLevel.ALL) {
      mDecomposer = new DoNothingDecomposer();
    } else if (trimSplitType == DecomposerType.TRIM || mParams.ionTorrent()) {
      // XXX I'm not sure why we do trimming only for ion torrent -- perhaps this could now use our ordinary trim/split
      mDecomposer = new Trimmer(variantAlleleTrigger);
    } else if (trimSplitType == DecomposerType.TRIMSPLIT) {
      mDecomposer = new SimpleDecomposer(mConfig.getDenovoChecker(), variantAlleleTrigger);
    } else {
      assert trimSplitType == DecomposerType.ALIGN;
      mDecomposer = new AligningDecomposer(mConfig.getDenovoChecker(), variantAlleleTrigger);
    }
    mAnnotators.addAll(mConfig.getVcfAnnotators());
    // AVR annotator comes last because it wants to use other annotations
    if (mParams.avrModelFile() != null) {
      Diagnostic.userLog("Loading AVR model: " + mParams.avrModelFile());
      mAnnotators.add(new ModelFactory(mParams.avrModelFile(), mParams.minAvrScore()).getModel());
    }
    mFilters.addAll(mConfig.getVcfFilters());

    String[] genomeNames = mConfig.getGenomeNames();
    if (genomeNames.length == 1) {
      genomeNames = new String[] {}; // Special case for singleton caller, map all records to 0
    }
    final SingletonPopulatorFactory<VariantAlignmentRecord> pf = new SingletonPopulatorFactory<>(new VariantAlignmentRecordPopulator(MultisampleUtils.chooser(mParams), mParams.minBaseQuality(), genomeNames));
    mWrapper = new ThreadedMultifileIteratorWrapper<>(new SamReadingContext(mParams.mapped(), mParams.ioThreads(), mParams.filterParams(), mParams.uberHeader(), mReferenceSequences), pf);
    final SAMSequenceDictionary dict = mWrapper.header().getSequenceDictionary();
    mSequences = dict.getSequences();
    if (mSequences.size() < 1) {
      throw new NoTalkbackSlimException("SAM file does not contain sequence dictionary.");
    }
    mFormatter = mConfig.getOutputFormatter(mParams);
    mVcfHeader = mFormatter.makeHeader(mParams, mParams.uberHeader());
    for (final VcfAnnotator annot : mAnnotators) {
      annot.updateHeader(mVcfHeader);
    }
    for (final VcfFilter filter : mFilters) {
      filter.setHeader(mVcfHeader);
    }
    mOut = new VcfWriterFactory().async(true).zip(mParams.blockCompressed()).make(mVcfHeader, mParams.vcfFile());
    mOut = new StatisticsVcfWriter<>(mOut, mStatistics);
    mBedFilterRegions = (mParams.regionsFilterBedFile() == null) ? null : BedUtils.regions(mParams.regionsFilterBedFile());

    Diagnostic.developerLog("Chunk size is " + mParams.chunkSize());
    Diagnostic.developerLog("Lookahead is " + mParams.threadingLookAhead());
  }

  @Override
  @SuppressWarnings("try")
  protected void exec() throws IOException {
    try {
      SamUtils.checkUberHeaderAgainstReference(mReferenceSequences, mParams.uberHeader(), !mParams.ignoreIncompatibleSamHeaders());
      init();
      final Map<String, Long> sequenceNameMap = ReaderUtils.getSequenceNameMap(mReferenceSequences);
      for (final SAMSequenceRecord r : mSequences) {
        final String sequenceName = r.getSequenceName(); //mReferenceSequences.name(l);
        // Only process this sequence if we are doing them all, or if it is
        // in the restriction specified by the user
        if (!mWrapper.context().hasRegions() || mWrapper.context().referenceRanges().containsSequence(sequenceName)) {
          mWrapper.setSequenceId(r.getSequenceIndex());
          if (mWrapper.hasNext()) {
            if (!sequenceNameMap.containsKey(sequenceName)) { //this means our SDF does not have reference
              throw new NoTalkbackSlimException("Reference SDF does not contain sequence '" + sequenceName + "'");
            }
            final long sdfSeqId = sequenceNameMap.get(sequenceName);
            final int sequenceLength = mReferenceSequences.length(sdfSeqId);
            final byte[] sequenceNt = new byte[sequenceLength];
            mReferenceSequences.read(sdfSeqId, sequenceNt);
            processAnEntireSequence(sequenceName, sequenceNt);
          }
        }
      }
      logRecordCounts();
      mStatistics.setExcessiveCoverageCount(mExcessiveCoverageCount);
      mStatistics.setExcessiveHypothesesCount(mExcessiveHypothesesCount);
      mStatistics.setNoHypothesesCount(mNoHypothesesCount);

      final long totalCalls = mStatistics.getTotalFiltered() + mStatistics.getTotalPassed();
      final double excessCoverageFraction = 100.0 * mExcessiveCoverageCount / totalCalls;
      if (mExcessiveCoverageCount > MIN_CALLS_FOR_COVERAGE_WARNING && excessCoverageFraction > COVERAGE_WARNING_THRESHOLD) {
        Diagnostic.warning("A large fraction of sites had coverage much higher than expected!  Check that input alignments have been calibrated with correct regions or that an appropriate --" + AbstractMultisampleCli.COVERAGE_BYPASS_FLAG + " value is supplied.");
      }

    } finally {
      try (VcfWriter ignored1 = mOut;
           OutputStream ignored2 = mBedOut;
           ThreadedMultifileIteratorWrapper<VariantAlignmentRecord> ignored3 = mWrapper;
           SequencesReader ignored4 = mReferenceSequences) {  // We want the sexy closing side effects.
        if (mConfig != null) {
          mConfig.close();
        }
      }
    }
    if (mParams.blockCompressed() && mParams.outputIndex()) {
      VcfUtils.createVcfTabixIndex(mParams.vcfFile());
      BedUtils.createBedTabixIndex(mParams.bedFile());
    }
  }

  @SuppressWarnings("unchecked")
  private static  List<VcfRecord> getVcfList(final Object o) {
    return (List<VcfRecord>) o;
  }
}
