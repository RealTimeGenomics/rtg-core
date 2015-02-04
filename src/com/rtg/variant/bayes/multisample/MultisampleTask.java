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
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import com.rtg.bed.BedUtils;
import com.rtg.launcher.ParamsTask;
import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.reference.Ploidy;
import com.rtg.reference.SexMemo;
import com.rtg.relation.Family;
import com.rtg.relation.Relationship;
import com.rtg.relation.Relationship.RelationshipType;
import com.rtg.sam.CircularBufferMultifileSinglePassReaderWindowSync;
import com.rtg.sam.ReaderRecord;
import com.rtg.sam.ReaderWindow;
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
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.VariantStatistics;
import com.rtg.variant.avr.ModelFactory;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.complex.IonTorrentCallFilter;
import com.rtg.variant.bayes.complex.Trimming;
import com.rtg.variant.bayes.multisample.cancer.SomaticCallerConfiguration;
import com.rtg.variant.bayes.multisample.family.DiseasedFamilyCallerConfiguration;
import com.rtg.variant.bayes.multisample.family.FamilyCallerConfiguration;
import com.rtg.variant.bayes.multisample.multithread.DependenciesMultiSample;
import com.rtg.variant.bayes.multisample.multithread.EventListMultiSample;
import com.rtg.variant.bayes.multisample.multithread.JobIdMultisample;
import com.rtg.variant.bayes.multisample.multithread.MultisampleStatistics;
import com.rtg.variant.bayes.multisample.population.PopulationCallerConfiguration;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.format.VariantOutputVcfFormatter;
import com.rtg.variant.util.VariantUtils;
import com.rtg.vcf.VcfAnnotator;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;
import com.rtg.vcf.VcfWriter;
import com.rtg.vcf.header.VcfHeader;

import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

/**
 */
public class MultisampleTask extends ParamsTask<VariantParams, VariantStatistics> {

  private static final int ION_TORRENT_HYPER_COMPLEX = 21;

  private final JobStatistics<JobIdMultisample> mJobStatistics = new MultisampleStatistics();
  private final OutputStream mBedOut;
  private final SexMemo mSexMemo;
  private final List<VcfAnnotator> mAnnotators = new ArrayList<>();
  protected final SequencesReader mReferenceSequences;
  private AbstractJointCallerConfiguration mConfig;
  private VcfWriter mOut;
  private VariantOutputVcfFormatter mFormatter;
  private VcfHeader mVcfHeader;
  private ThreadedMultifileIteratorWrapper<VariantAlignmentRecord> mWrapper;
  private List<SAMSequenceRecord> mSequences;
  private ReferenceRegions mBedFilterRegions;
  private ParallelProgress mPP = null;

  private final JointCallerConfigurator mConfigurator;

  // This configurator tries to guess the appropriate configuration based on the structure of the pedigree
  private static final JointCallerConfigurator DEFAULT_CONFIGURATOR = new JointCallerConfigurator() {
    @Override
    public AbstractJointCallerConfiguration getConfig(VariantParams params, String[] outputSampleNames) throws IOException {
      final Relationship[] derived = params.genomeRelationships().relationships(RelationshipType.ORIGINAL_DERIVED);
      if (derived.length > 0) {
        final int numberOfGenomes = params.genomeRelationships().genomes().length;
        if (derived.length != 1 || numberOfGenomes != 2) {
          throw new RuntimeException("Single original-derived relationship required in genome relationship file.");
        }
        return new SomaticCallerConfiguration.Configurator().getConfig(params, outputSampleNames);
      } else {
        final Relationship[] parentChild = params.genomeRelationships().relationships(RelationshipType.PARENT_CHILD);
        if (parentChild.length > 0) {
          // We're in family calling mode
          final Family family = Family.getFamily(params.genomeRelationships());
          if (family.isOneParentDiseased()) {
            return new DiseasedFamilyCallerConfiguration.Configurator().getConfig(params, outputSampleNames);
          } else {
            return new FamilyCallerConfiguration.Configurator().getConfig(params, outputSampleNames);
          }
        } else {
          // We're in population calling mode
          return new PopulationCallerConfiguration.Configurator().getConfig(params, outputSampleNames);
        }
      }
    }
  };

  /**
   * @param params command line parameters.
   * @param configurator the joint caller configuration creator
   * @param defaultOutput where the calling results are written.
   * @param statistics the variant statistics collector
   * @param usageMetric count of the number of records read.
   * @throws IOException when writing output.
   */
  public MultisampleTask(VariantParams params, JointCallerConfigurator configurator, OutputStream defaultOutput, VariantStatistics statistics, UsageMetric usageMetric) throws IOException {
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

  /**
   * @param params command line parameters.
   * @param defaultOutput where the calling results are written.
   * @param statistics the variant statistics collector
   * @param usageMetric to be updated with the count of the number of nucleotides read.
   * @throws IOException when writing output.
   */
  public MultisampleTask(VariantParams params, OutputStream defaultOutput, VariantStatistics statistics, UsageMetric usageMetric) throws IOException {
    this(params, DEFAULT_CONFIGURATOR, defaultOutput, statistics, usageMetric);
  }

  private List<ModelInterface<?>> stepModels(final IndividualSampleProcessor<?>[] ssProcessors, int pos) {
    final List<ModelInterface<?>> models = new ArrayList<>();
    for (final IndividualSampleProcessor<?> ssProcessor : ssProcessors) {
      models.add(ssProcessor.step(pos));
    }
    return models;
  }
  private void mindlesslyAdvanceIndels(final IndividualSampleProcessor<?>[] ssProcessors, int pos) {
    for (final IndividualSampleProcessor<?> ssProcessor : ssProcessors) {
      ssProcessor.stepIndel(pos);
    }
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

        for (int pos = start; pos < end; ) {
          if (statusInterval.contains(pos)) {
            final byte status = statusInterval.get(pos);
            final int oldpos = pos;
            do {
              stepModels(ssProcessors, pos);
              mindlesslyAdvanceIndels(ssProcessors, pos);
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
            final List<ModelInterface<?>> models = stepModels(ssProcessors, pos);

            int refHyp = Integer.MIN_VALUE;
            for (ModelInterface<?> m : models) {
              if (m.hypotheses().size() > 0) {
                refHyp = m.reference(); /// XXX This is still a bit icky. We are assuming that if you have hypotheses all references are the same (maybe not true in cancer)
              }
            }
            if (refHyp != Integer.MIN_VALUE) {
              final HaploidDiploidHypotheses<HypothesesPrior<Description>> hypotheses = mConfig.getSnpHypotheses(refHyp, refName, pos);
              final Variant variant = jointCaller.makeCall(refName, pos, pos + 1, template, models, hypotheses);
              if (variant != null) {
                calls.add(variant);
              }
            }
            pos++;
          }
        }
      }
    }
    //TODO move to the output task so come out in sequence
    mPP.updateProgress(chunkInfo.percent(end));
    return maxReadLen;
  }

  // Sets the status of any positions within the interval to SKIP if they are contained within a no-call range entry (one without metadata)
  private static void addRangeStatuses(StatusInterval statusInterval, List<RangeList.RangeData<String>> ranges, int startIndex, int endPos) {
    for (int rangeIndex = startIndex; rangeIndex < ranges.size(); rangeIndex++) {
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
    mInvalidRecords++;
    if (mInvalidRecords <= 5) {
      Diagnostic.warning(WarningType.SAM_BAD_FORMAT_WARNING1, record.toString());
      Diagnostic.userLog("Invalid record: " + record.toString());
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

    JobFactoryMultiSample(final ChunkInfo info, final String refName, final byte[] refNts) throws IOException {
      mInfo = info;
      String[] genomeNames = mConfig.getGenomeNames();
      if (genomeNames.length == 1) {
        genomeNames = new String[] {}; // Special case for singleton caller, map all records to 0
      }
      final int depth = mParams.maxCoverageBypass().thresholdTotal(refName);
      final VariantAlignmentRecordPopulator pop = new VariantAlignmentRecordPopulator(genomeNames);
      final RegionRestriction restriction = new RegionRestriction(refName, info.start(), info.end());

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
      public MultisampleJob(JobIdMultisample id, Result[] arguments) {
        super(id);
        mArguments = arguments;
      }

      @Override
      public String toString() {
        final int start = id().time() * mInfo.chunkSize() + mInfo.start();
        final int end = Math.min(start + mInfo.chunkSize(), mInfo.end());
        return super.toString() + " " + mRefName + ":" + start + "-" + end;
      }

      protected List<Variant> getList(Object o) {
        @SuppressWarnings("unchecked")
        final List<Variant> ret = (List<Variant>) o;
        return ret;
      }
    }

    private class IncrJob extends MultisampleJob {
      public IncrJob(JobIdMultisample id, Result[] arguments) {
        super(id, arguments);
      }

      @Override
      public Result run() throws IOException {
        final int start = id().time() * mInfo.chunkSize() + mInfo.start();
        final int end = Math.min(start + mInfo.chunkSize(), mInfo.end());
        final List<Variant> calls = new ArrayList<>();
        final int maxReadLen = processNtPositions(calls, mJointCaller, mInfo, mRefNts, mBuffer, start, end);
        final boolean simpleRepeats = mParams.simpleRepeatExtension() && !mParams.ionTorrent();
        final Complexities cx = new Complexities(calls, mRefName, start, end, mParams.interestingSeparation(), mParams.hyperComplexLength(), mRefNts, simpleRepeats, mConfig.getSiteSpecificPriors());
        return new Result(cx, maxReadLen);
      }
    }

    private class DanglingJob extends MultisampleJob {
      public DanglingJob(JobIdMultisample id, Result[] arguments) {
        super(id, arguments);
      }

      @Override
      public Result run() {
        final Complexities last = mArguments[0] == null ? null : (Complexities) mArguments[0].result(0);
        final Complexities cx = mArguments[1] == null ? null : (Complexities) mArguments[1].result(0);
        Complexities.fixDangling(last, cx);
        return new Result(last);
      }
    }

    private class BedJob extends MultisampleJob {
      public BedJob(JobIdMultisample id, Result[] arguments) {
        super(id, arguments);
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
      public ComplexJob(JobIdMultisample id, Result[] arguments) {
        super(id, arguments);
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
              if (!v.isIndel() && !v.isOverflow() && (v.getLocus().getStart() >= complexRegions.startOfChunk() && v.getLocus().getStart() < complexRegions.endOfChunk())) {
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
          mBuffer.flush(complexRegions.startOfChunk(), complexRegions.endOfChunk());
          return new Result(calls);
        } else {
          return new Result((Object) null);
        }
      }
    }

    private class FilterJob extends MultisampleJob {
      public FilterJob(JobIdMultisample id, Result[] arguments) {
        super(id, arguments);
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
          final EquivalentFilter filter = new EquivalentFilter(mRefNts, prevLastCall);
          equivFiltered = filter.filter(split, maxReadLen);
          lastCalls = filter.lastCall(); // Remember last calls for checking equivalence across chunks
        }

        final List<Variant> filtered = locusAndIonTorrentFilters(equivFiltered, mRefNts);
        return new Result(filtered, lastCalls);
      }
    }

    private class OutJob extends MultisampleJob {
      public OutJob(JobIdMultisample id, Result[] arguments) {
        super(id, arguments);
      }

      @Override
      public Result run() throws IOException {
        final List<Variant> filtered = getList(mArguments[1].result(0));
        if (filtered != null) {
          for (final Variant v : filtered) {
            final VcfRecord record = mFormatter.makeVcfRecord(v);
            for (VcfAnnotator annot : mAnnotators) {
              annot.annotate(record);
            }
            mStatistics.tallyVariant(mVcfHeader, record);
            mOut.write(record);
          }
        }
        if (mBuffer.finishedTo() < Math.min(id().time() * mInfo.chunkSize() + mInfo.start(), mInfo.end())) { //flushing should be keeping up with output
          throw new RuntimeException("Failed to flush chunk: " + mBuffer.finishedTo() + " : " + id().time() * mInfo.chunkSize());
        }
        return null;
      }
    }
  }

  private List<Variant> trimSplit(List<Variant> merged) {
    if (!mParams.enableTrimSplit()) {
      return merged;
    }
    final List<Variant> calls = new ArrayList<>();
    for (Variant variant : merged) {
      if (variant.hasCallNames()) {
        calls.addAll(Trimming.trimSplit(mParams, variant, mConfig.getDenovoCorrector()));
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
    assert ranges.size() > 0;
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
      final Scheduler<JobIdMultisample> sched = new SchedulerSynchronized<>(depen, jobFac, eventList, null, mJobStatistics, mParams.threadingLookAhead());
      //final Scheduler<JobIdMultisample> sched = new SchedulerSynchronized<>(depen, jobFac, eventList, System.err, mJobStatistics, mParams.threadingLookAhead());
      final String msg = "Processing " + refName;
      final Executor<JobIdMultisample> exec = createExecutor(sched, msg, mParams);
      exec.run();
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

    final String invalidRecordsWarning = invalidRecords + " records skipped because of SAM format problems.";
    if (invalidRecords > 0) {
      Diagnostic.warning(invalidRecordsWarning);
    } else {
      Diagnostic.userLog(invalidRecordsWarning);
    }
    final long filteredRecords = mWrapper.getFilteredRecordsCount();
    if (filteredRecords > 0) {
      Diagnostic.userLog(filteredRecords + " records skipped due to input filtering criteria");
    }
    if (mTossedRecords > 0) {
      Diagnostic.userLog(mTossedRecords + " records skipped in extreme coverage regions");
    }
    final long numDeduped = mWrapper.getDuplicateRecordsCount();
    if (numDeduped != 0) {
      Diagnostic.userLog(numDeduped + " records skipped due to duplicate detection");
    }
    Diagnostic.userLog(mWrapper.getOutputRecordsCount() + " records processed");
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
      Diagnostic.developerLog("New params:" + mParams.toString());
    }

    // Output for every sample for which mappings are provided
    final Set<String> samples = new TreeSet<>();
    for (final SAMReadGroupRecord ih : mParams.uberHeader().getReadGroups()) {
      final String sample = ih.getSample();
      if (sample != null) {
        samples.add(sample);
      }
    }

    mConfig = mConfigurator.getConfig(mParams, samples.toArray(new String[samples.size()]));
    mAnnotators.addAll(mConfig.getVcfAnnotators());

    if (mParams.avrModelFile() != null) {
      Diagnostic.userLog("Loading AVR model: " + mParams.avrModelFile());
      mAnnotators.add(new ModelFactory(mParams.avrModelFile(), mParams.minAvrScore()).getModel());
    }

    String[] genomeNames = mConfig.getGenomeNames();
    if (genomeNames.length == 1) {
      genomeNames = new String[] {}; // Special case for singleton caller, map all records to 0
    }
    final SingletonPopulatorFactory<VariantAlignmentRecord> pf = new SingletonPopulatorFactory<>(new VariantAlignmentRecordPopulator(genomeNames));
    mWrapper = new ThreadedMultifileIteratorWrapper<>(mParams.mapped(), mParams.ioThreads(), pf, mParams.filterParams(), mParams.uberHeader());
    final SAMSequenceDictionary dict = mWrapper.header().getSequenceDictionary();
    mSequences = dict.getSequences();
    if (mSequences.size() < 1) {
      throw new NoTalkbackSlimException("SAM file does not contain sequence dictionary.");
    }
    mFormatter = mConfig.getOutputFormatter(mParams);
    mVcfHeader = mFormatter.makeHeader(mParams, mParams.uberHeader());
    for (VcfAnnotator annot : mAnnotators) {
      annot.updateHeader(mVcfHeader);
    }
    mOut = new VcfWriter(mVcfHeader, mParams.vcfStream());
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
          if (!sequenceNameMap.containsKey(sequenceName)) { //this means our SDF does not have reference
            throw new NoTalkbackSlimException("Reference SDF does not contain sequence '" + sequenceName + "'");
          }
          final long sdfSeqId = sequenceNameMap.get(sequenceName);
          final int sequenceLength = mReferenceSequences.length(sdfSeqId);
          mWrapper.setSequenceId(r.getSequenceIndex());
          if (mWrapper.hasNext()) {
            assert sequenceLength == r.getSequenceLength() : "" + sequenceLength + " != " + r.getSequenceLength();
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
}
