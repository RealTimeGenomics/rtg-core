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

package com.rtg.metagenomics;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.launcher.ParamsTask;
import com.rtg.launcher.ReaderParams;
import com.rtg.metagenomics.krona.KronaSpeciesNode;
import com.rtg.metagenomics.krona.KronaSpeciesReportWriter;
import com.rtg.metagenomics.matrix.Vector;
import com.rtg.reader.PrereadNamesInterface;
import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.sam.SamReadingContext;
import com.rtg.sam.SamRecordPopulator;
import com.rtg.sam.SamUtils;
import com.rtg.sam.ThreadedMultifileIterator;
import com.rtg.taxonomy.TaxonNode;
import com.rtg.taxonomy.Taxonomy;
import com.rtg.taxonomy.TaxonomyUtils;
import com.rtg.usage.UsageMetric;
import com.rtg.util.BoundedDouble;
import com.rtg.util.ComparablePair;
import com.rtg.util.Environment;
import com.rtg.util.HtmlReportHelper;
import com.rtg.util.IORunnable;
import com.rtg.util.ReverseDouble;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.SingletonPopulatorFactory;
import com.rtg.util.Utils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.LineWriter;

import htsjdk.samtools.SAMRecord;

/**
 */
@TestClass("com.rtg.metagenomics.SpeciesCliTest")
class SpeciesTask extends ParamsTask<SpeciesParams, SpeciesStatistics> {

  private static final String TAB = "\t";
  private static final String COMMENT_CHAR = "#";

  /** The attribute used in a comment line to indicate the SDF ID of the template */
  private static final String TEMPLATE_SDF_ATTRIBUTE = "TEMPLATE-SDF-ID";
  private static final int MAX_WARNING = 5;


  private static final String ABUNDANCE_LABEL = "abundance";
  private static final String ABUNDANCE_LOW_LABEL = "abundance-low";
  private static final String ABUNDANCE_HIGH_LABEL = "abundance-high";
  private static final String DNA_LABEL = "DNA-fraction";
  private static final String DNA_LOW_LABEL = "DNA-fraction-low";
  private static final String DNA_HIGH_LABEL = "DNA-fraction-high";
  private static final String CONFIDENCE_LABEL = "confidence";
  private static final String COVERAGE_DEPTH_LABEL = "coverage-depth";
  private static final String COVERAGE_BREADTH_LABEL = "coverage-breadth";
  private static final String GENOME_LENGTH_LABEL = "reference-length";
  private static final String MAPPED_READS_LABEL = "mapped-reads";
  private static final String HAS_SEQUENCE_LABEL = "has-reference";
  private static final String SEQUENCE_COUNT_LABEL = "taxa-count";
  private static final String TAXON_ID_LABEL = "taxon-id";
  private static final String PARENT_ID_LABEL = "parent-id";
  private static final String RANK_LABEL = "rank";
  private static final String TAXONOMY_NAME_LABEL = "taxonomy-name";

  private static final String VERSION = "Version" + TAB + Environment.getVersion() + ", species v2.1";
  static final String SPECIES_HEADER = ABUNDANCE_LABEL + TAB + ABUNDANCE_LOW_LABEL + TAB + ABUNDANCE_HIGH_LABEL
      + TAB + DNA_LABEL + TAB + DNA_LOW_LABEL + TAB + DNA_HIGH_LABEL
      + TAB + CONFIDENCE_LABEL + TAB + COVERAGE_DEPTH_LABEL + TAB + COVERAGE_BREADTH_LABEL + TAB + GENOME_LENGTH_LABEL
      + TAB + MAPPED_READS_LABEL + TAB + HAS_SEQUENCE_LABEL + TAB + SEQUENCE_COUNT_LABEL + TAB + TAXON_ID_LABEL + TAB + PARENT_ID_LABEL
      + TAB + RANK_LABEL + TAB + TAXONOMY_NAME_LABEL;

  private static final String NUM_FORMAT = "%1.4g";
  private static final String NUM_FORMAT_CONFIDENCE = "%1.2g";

  protected final Map<String, Integer> mSequenceMap = new HashMap<>();
  private final Map<String, ArrayList<Integer>> mHits = new HashMap<>();
  private final Set<Integer> mSpeciesWithHits = new HashSet<>();

  //TODO get rid of these protected declarations - this is not a civilized way to test things
  //or to initialize them. Maybe they should be pulled out into a separate class
  protected double mSampleFactor = 1.0;
  final Map<Integer, SpeciesInfo> mSpeciesInfo = new HashMap<>();
  protected SpeciesMap mSpeciesMap = null;
  protected Taxonomy mTaxonomy;

  static class MyInteger {
    int mValue = 1;
  }

  public SpeciesTask(final SpeciesParams params, final OutputStream out, final UsageMetric usageMetric) {
    super(params, out, new SpeciesStatistics(params.directory()), usageMetric);
  }

  static SpeciesMap getSpeciesMap(Map<String, Integer> sequenceMap, SequencesReader sr, Taxonomy taxonomy) throws IOException {
    final SpeciesMap speciesMap = new SpeciesMap();
    // Preload species map to ensure names in correct order
    final PrereadNamesInterface names = sr.names();
    for (int k = 0; k < names.length(); k++) {
      final String name = names.name(k);
      final Integer taxonId = sequenceMap.get(name);
      TaxonNode taxon = taxonomy.get(taxonId);
      do {
        speciesMap.id(taxon.getId());
        taxon = taxon.getParent();
      } while (taxon != null && taxon.getId() != Taxonomy.ROOT_ID);
    }

    return speciesMap;
  }

  protected void buildDefaultTaxonomy(File referenceMap, SequencesReader sr) throws IOException {
    mTaxonomy = new Taxonomy();
    int warningCount = 0;
    mTaxonomy.addNode(Taxonomy.ROOT_ID, -1, "root", "no rank");

    int nextId = 2;
    final Map<String, Long> namemap = ReaderUtils.getSequenceNameMap(sr);
    if (referenceMap != null) {
      try (BufferedReader inputList = new BufferedReader(new FileReader(referenceMap))) {
        String line;
        while ((line = inputList.readLine()) != null) {
          if (line.charAt(0) == SpeciesCli.COMMENT_CHAR) {
            continue;
          }
          final int tabIndex = line.indexOf('\t');
          if (tabIndex == -1 || tabIndex == 0 || tabIndex == line.length() - 1) {
            throw new NoTalkbackSlimException("The input file " + referenceMap.getPath() + " for relabel-species-file flag contains invalid entry. line: " + line);
          }
          final String seqName = line.substring(0, tabIndex).trim();
          final String speciesName = line.substring(tabIndex).trim();
          if (!namemap.containsKey(seqName)) {
            if (warningCount++ < MAX_WARNING) {
              Diagnostic.warning("Template SDF does not contain any sequence named \"" + seqName + "\" specified in the " + referenceMap.getPath() + " file");
            }
            if (warningCount == MAX_WARNING) {
              Diagnostic.warning("Subsequent warnings of this type will not be shown.");
            }
          }
          Integer taxonId = mSequenceMap.get(speciesName);
          if (taxonId == null) {
            taxonId = nextId++;
            mSequenceMap.put(speciesName, taxonId);
            mTaxonomy.addNode(taxonId, 1, speciesName, "species");
          }
          mSequenceMap.put(seqName, taxonId);
        }
      }
    }
    final int taxonIdBefore = nextId;
    for (long i = 0; i < sr.numberSequences(); i++) {
      final String shortName = sr.name(i);
      if (!mSequenceMap.containsKey(shortName)) {
        final Integer taxonId = nextId++;
        mSequenceMap.put(shortName, taxonId);
        mTaxonomy.addNode(taxonId, 1, sr.fullName(i), "species");
      }
    }
    if (nextId > taxonIdBefore) {
      Diagnostic.userLog("There were " + (nextId - taxonIdBefore) + " reference sequences not relabeled");
    }
    if (warningCount >= MAX_WARNING) {
      Diagnostic.warning("There were " + warningCount + " names not present in the template SDF.");
    }
  }

  static HashMap<String, String> getRenameMap(File referenceMap, SequencesReader sr) throws IOException {
    final Map<String, Long> namemap = ReaderUtils.getSequenceNameMap(sr);
    int warningCount = 0;
    final HashMap<String, String> renameMap = new HashMap<>();
    if (referenceMap != null) {
      try (BufferedReader inputList = new BufferedReader(new FileReader(referenceMap))) {
        String line;
        while ((line = inputList.readLine()) != null) {
          if (line.charAt(0) == SpeciesCli.COMMENT_CHAR) {
            continue;
          }
          final int tabIndex = line.indexOf('\t');
          if (tabIndex == -1 || tabIndex == 0 || tabIndex == line.length() - 1) {
            throw new NoTalkbackSlimException("The input file " + referenceMap.getPath() + " for relabel-species-file flag contains invalid entry. line: " + line);
          }
          final String referenceName = line.substring(0, tabIndex).trim();
          final String speciesName = line.substring(tabIndex).trim();
          if (!namemap.containsKey(referenceName)) {
            if (warningCount++ < MAX_WARNING) {
              Diagnostic.warning("Template SDF does not contain any sequence named \"" + referenceName + "\" specified in the " + referenceMap.getPath() + " file");
            }
            if (warningCount == MAX_WARNING) {
              Diagnostic.warning("Subsequent warnings of this type will not be shown.");
            }
          }
          renameMap.put(referenceName, speciesName);
        }
      }
    }
    for (long i = 0; i < sr.numberSequences(); i++) {
      final String shortName = sr.name(i);
      final String suffix = sr.nameSuffix(i);
      if (!renameMap.containsKey(shortName) && suffix.length() != 0) {
        renameMap.put(shortName, sr.fullName(i));
      }
    }
    if (warningCount >= MAX_WARNING) {
      Diagnostic.warning("There were " + warningCount + " names not present in the template SDF.");
    }
    return renameMap;
  }

  @Override
  protected void exec() throws IOException {
    Diagnostic.progress("Reading names.");
    final SequencesReader sr = mParams.genome().reader();

    if (TaxonomyUtils.hasTaxonomyInfo(sr)) {
      mSequenceMap.putAll(TaxonomyUtils.loadTaxonomyMapping(sr));
      final Taxonomy baseTaxonomy = TaxonomyUtils.loadTaxonomy(sr);
      if (!baseTaxonomy.isConsistent()) {
        throw new NoTalkbackSlimException("The taxonomy in the provided SDF is invalid: " + baseTaxonomy.getInconsistencyReason());
      }
      mTaxonomy = baseTaxonomy.subset(mSequenceMap.values());

      // ensure all sequences accounted for
      final PrereadNamesInterface names = sr.names();
      for (int k = 0; k < names.length(); k++) {
        final String name = names.name(k);
        if (!mSequenceMap.containsKey(name)) {
          throw new NoTalkbackSlimException("Reference SDF contains sequences not referenced by sequence lookup");
        }
      }
    } else {
      buildDefaultTaxonomy(mParams.referenceMap(), sr);
    }
    Diagnostic.userLog("Number of reference sequences: " + mSequenceMap.size());
    Diagnostic.userLog("Number of reference species: " + new HashSet<>(mSequenceMap.values()).size());
    Diagnostic.userLog("Number of taxonomy nodes: " + mTaxonomy.size());
    mSpeciesMap = getSpeciesMap(mSequenceMap, sr, mTaxonomy);

    SamUtils.checkUberHeaderAgainstReference(sr, SamUtils.getUberHeader(mParams.mapped(), true, null), false);

    accumulateMappings(sr);

  }

  private void processMappings(final SequencesReader sr) throws IOException {
    final Frag[] frags = preprocessHits();

    Diagnostic.progress("Block Construction Started");
    final BlockMapping blockMapping = new BlockMapping(frags, mSpeciesMap.taxonIds().length);
    Diagnostic.developerLog(blockMapping.statistics());

    final long[] genomeLengths = new long[mSpeciesMap.taxonIds().length];

    final ReaderParams genomes = mParams.genome();
    final int[] initGenomeLengths = genomes.lengths();
    for (int i = 0; i < initGenomeLengths.length; i++) {
      final int id = mSpeciesMap.id(mSequenceMap.get(genomes.reader().name(i)));
      genomeLengths[id] += initGenomeLengths[i];
    }
    //separate into blocks
    final BlockInfo blockInfo = new BlockInfo(-1, null, frags, mSpeciesMap, genomeLengths, mParams.verbose());
    Diagnostic.progress("Block Construction Finished");
    final BlockResult result = solveBlocks(frags, blockMapping, genomeLengths, blockInfo);
    //OUTPUT
    try (final LineWriter out = new LineWriter(new OutputStreamWriter(mParams.speciesStream()))) {
      result(sr, out, result, blockInfo);
    }
  }

  private Frag[] preprocessHits() {
    Diagnostic.progress("Pre-processing Started");
    final Frag[] frags;
    // Collapse identical frags
    final HashMap<Frag, SpeciesTask.MyInteger> uniq = new HashMap<>();
    for (final ArrayList<Integer> g : mHits.values()) {
      Collections.sort(g); // needed for frag construction?
      final Frag f = new Frag(g);
      if (uniq.containsKey(f)) {
        uniq.get(f).mValue++;
      } else {
        uniq.put(f, new SpeciesTask.MyInteger());
      }
    }
    // Zero-frequency adjustment. Add one fragment per species. Doing this helps with
    // convergence and confidence evaluation later on.
    for (final Integer taxonId : new HashSet<>(mSequenceMap.values())) {
      final int species = mSpeciesMap.id(taxonId);
      if (!mSpeciesWithHits.contains(species)) {
        continue;
      }
      //System.err.println("taxonId=" + taxonId + " species=" + species);
      final Frag f = new Frag(Collections.singletonList(species));
      if (uniq.containsKey(f)) {
        uniq.get(f).mValue++;
      } else {
        uniq.put(f, new SpeciesTask.MyInteger());
      }
    }
    frags = new Frag[uniq.size()];
    int k = 0;
    for (final Map.Entry<Frag, SpeciesTask.MyInteger> e : uniq.entrySet()) {
      final Frag f = e.getKey();
      f.setMultiplicity(e.getValue().mValue);
      frags[k++] = f;
    }
    Diagnostic.developerLog("Frag collapsed count = " + frags.length + ", original count = " + mHits.size());
    Diagnostic.progress("Pre-processing Finished");
    return frags;
  }

  protected void accumulateMappings(final SequencesReader sr) throws IOException {
    Diagnostic.progress("SAM Reading Started");
    String lastSequenceName = null;
    double cov = 0.0;
    long bcov = 0;
    long bcount = 0;
    long oldPos = 0;
    long usageStats = 0L;
    double mappedSpec = 0.0;
    double mappedReads = 0.0;
    double unmappedReads = 0.0;
    final SamReadingContext context = new SamReadingContext(mParams.mapped(), 1, mParams.filterParams(), SamUtils.getUberHeader(mParams.mapped()));
    try (final ThreadedMultifileIterator<SAMRecord> it = new ThreadedMultifileIterator<>(context, new SingletonPopulatorFactory<>(new SamRecordPopulator()))) {
      while (it.hasNext()) {
        usageStats++;
        final SAMRecord rec = it.next();
        final Integer ih = rec.getIntegerAttribute(SamUtils.ATTRIBUTE_IH);
        final double mappedIncr = (ih != null && ih > 0) ? 1.0 / ih : 1.0;
        if (rec.getReadUnmappedFlag()) {
          unmappedReads += 1.0;
        } else {
          mappedReads += mappedIncr;
          final String readId = rec.getReadName();
          final String sequenceName = rec.getReferenceName();
          final Integer taxonId = mSequenceMap.get(sequenceName);
          if (taxonId == null) {
            // something wrong - maybe mappings against different reference ?
            Diagnostic.warning("Could not find taxon ID for sequence: " + sequenceName);
          } else {
            if (!sequenceName.equals(lastSequenceName)) {
              if (lastSequenceName != null) {
                final Integer lastTaxonId = mSequenceMap.get(lastSequenceName);
                mSpeciesInfo.put(lastTaxonId, new SpeciesInfo(cov, bcov + bcount, mappedSpec));
                if (mSpeciesInfo.containsKey(taxonId)) {
                  final SpeciesInfo info = mSpeciesInfo.get(taxonId);
                  cov = info.mCoverageDepth;
                  bcov = info.mCoverageBreadth;
                  mappedSpec = info.mMappedReads;
                } else {
                  cov = 0.0;
                  bcov = 0;
                  mappedSpec = 0.0;
                }
                oldPos = 0;
                bcount = 0;
              }
              lastSequenceName = sequenceName;
              Diagnostic.developerLog("Starting: " + sequenceName);
            }
            final Integer speciesId = mSpeciesMap.id(taxonId);
            ArrayList<Integer> hitList = mHits.get(readId);
            if (hitList == null) {
              hitList = new ArrayList<>();
              mHits.put(readId, hitList);
            }
            hitList.add(speciesId);
            mSpeciesWithHits.add(speciesId);
            final int len = rec.getReadLength();
            cov += len * mappedIncr;
            final int pos = rec.getAlignmentStart();
            if (pos != oldPos) {
              // Update breadth of coverage numbers
              final long delta = pos - oldPos;
              if (delta > bcount) {
                bcov += bcount;
                bcount = 0;
              } else {
                bcov += delta;
                bcount -= delta;
              }
              oldPos = pos;
            }
            bcount = Math.max(bcount, len);
            mappedSpec += mappedIncr;
          }
        }
      }
      if (lastSequenceName != null) {
        final Integer lastTaxonId = mSequenceMap.get(lastSequenceName);
        mSpeciesInfo.put(lastTaxonId, new SpeciesInfo(cov, bcov + bcount, mappedSpec));
      }
    }
    Diagnostic.userLog("Mapped reads: " + mappedReads);
    Diagnostic.userLog("Unmapped reads: " + unmappedReads);
    if (unmappedReads > 0) {
      mSampleFactor = mappedReads / (mappedReads + unmappedReads);
      Diagnostic.userLog("Fraction of sample mapped: " + mSampleFactor);
    } else {
      mSampleFactor = 1.0;
    }
    mUsageMetric.setMetric(usageStats);
    Diagnostic.progress("SAM Reading Finished");
    processMappings(sr);
  }

  private SpeciesResult solveBlocks(Frag[] frags, final BlockMapping blockMapping, long[] genomeLengths, final BlockInfo blockInfo) throws IOException {
    final BlockInfo[] subBlocks = blockMapping.subBlock(mParams, frags, mSpeciesMap, genomeLengths);
    final SubBlockResult[] subResults = new SubBlockResult[subBlocks.length];
    final SimpleThreadPool stp = new SimpleThreadPool(mParams.execThreads(), "SolveBlocks", true); // Block-level parallelism
    stp.enableBasicProgress(subBlocks.length);
    final ExecutorService pvalueExecutor = Executors.newFixedThreadPool(mParams.execThreads()); // Gives parallism to p-value calculation
    try {
      for (final BlockInfo subBlockInfo : subBlocks) {
        final IORunnable run = new IORunnable() {
          final BlockMapping mBlockMapping = blockMapping;
          final BlockInfo mCurrentBlock = subBlockInfo;
          final int mCurrentBlockId = subBlockInfo.id();
          final int[][] mMembersOf = new int[mSpeciesMap.size()][];
          final boolean[] mInBlock = new boolean[mSpeciesMap.size()];
          final Map<IdSet, Double> mKnownSolutions = Collections.synchronizedMap(new HashMap<IdSet, Double>());

          @Override
          public void run() throws IOException {
            //Diagnostic.progress("Starting block " + (mCurrentBlockId + 1) + "/" + subBlocks.length);
            Diagnostic.userLog("Starting block " + (mCurrentBlockId + 1));
            Diagnostic.developerLog(mCurrentBlock.toString());
            Diagnostic.developerLog(mCurrentBlock.fragInfo());

            // Set up block membership info
            for (int j = 0; j < mMembersOf.length; j++) {
              final TaxonNode taxon = mTaxonomy.get(mSpeciesMap.taxonId(j));
              final List<TaxonNode> members = taxon.depthFirstTraversal();
              mMembersOf[j] = new int[members.size()];
              int ji = 0;
              for (final TaxonNode t : members) {
                final Integer globalId = mSpeciesMap.get(t.getId());
                if (mBlockMapping.getBlockId(globalId) == mCurrentBlockId) {
                  mMembersOf[j][ji++] = mBlockMapping.getLocalId(globalId);
                }
              }
              mInBlock[j] = ji > 0;
              mMembersOf[j] = Arrays.copyOf(mMembersOf[j], ji);
            }

            subResults[mCurrentBlockId] = runBlock();

            for (int j = 0; j < mMembersOf.length; j++) {
              mMembersOf[j] = null;
            }

            Diagnostic.developerLog("Finished block " + mCurrentBlockId);
          }

          private SubBlockResult runBlock() throws IOException {
            final Species sp = new Species(mMembersOf, mCurrentBlock);
            final SubBlockResult subBlockResults = sp.solve(mCurrentBlock.getN() * mParams.minIter());

            // Calculate P values, for every global taxon id that could be affected:
            final Vector initialR = subBlockResults.getR();
            final List<Future<?>> pFutures = new ArrayList<>();
            for (int j = 0; j < mInBlock.length; j++) {
              if (mInBlock[j]) {
                final int nodeId = j;
                final Runnable run = new Runnable() {
                  @Override
                  public void run() {
                    final Species psp = new Species(mMembersOf, mCurrentBlock); // Species is not thread-safe, so make a new one
                    final IdSet idSet = new IdSet(mMembersOf[nodeId]);
                    final Double v = mKnownSolutions.get(idSet);
                    if (v == null) {
                      final double l = subBlockResults.getL();
                      final String name = mTaxonomy.get(mSpeciesMap.taxonId(nodeId)).getName();
                      Diagnostic.developerLog("B:" + mCurrentBlockId + " Doing p-value estimation for node " + nodeId + " " + name + " base L= " + l);
                      final double newL = psp.solveFixed(initialR, mMembersOf[nodeId], 2, l);
                      final double ll = newL - l;
                      Diagnostic.developerLog("B:" + mCurrentBlockId + " L0:" + l + " L':" + newL);
                      subBlockResults.getLikelihoods().set(nodeId, ll);
                      mKnownSolutions.put(idSet, ll);
                    } else {
                      subBlockResults.getLikelihoods().set(nodeId, v);
                    }
                  }
                };
                pFutures.add(pvalueExecutor.submit(run));
              }
            }
            // Now wait for all the p-values to be finished.
            for (final Future<?> f : pFutures) {
              try {
                f.get();
              } catch (final ExecutionException ie) {
                if (ie.getCause() instanceof NoTalkbackSlimException) {
                  throw (NoTalkbackSlimException) ie.getCause();
                } else if (ie.getCause() instanceof IOException) {
                  throw (IOException) ie.getCause();
                } else if (ie.getCause() instanceof Error) {
                  throw (Error) ie.getCause();
                }
                throw new RuntimeException(ie.getCause());
              } catch (final InterruptedException ie) {
                throw new NoTalkbackSlimException("Interrupted while calculating p-values.");
              }
            }
            return subBlockResults;
          }
        };
        stp.execute(run);
      }
      stp.terminate();
    } finally {
      pvalueExecutor.shutdownNow(); // Should already be finished by the time stp is finished.
    }

    Diagnostic.progress("Merging Block Results Started");
    Diagnostic.userLog("Merging block results.");
    final SpeciesResult result = new SpeciesResult(blockInfo.getN());
    for (int currentBlock = 0; currentBlock < subBlocks.length; currentBlock++) {
      blockMapping.mergeResults(result, subResults[currentBlock], currentBlock);
    }
    Diagnostic.progress("Merging Block Results Finished");
    return result;
  }

  private static double[] stdDevBounds(final int taxonId, final double wt, double v, double r) {
    // Expands out the estimate by 3 standard deviations, thus would expected to
    // encompass the correct result 99.7% of the time
    final double low; //= wt / es;
    final double high; //= wt * es;
    final double es = Math.exp(3 * Math.sqrt(v) / r);
    if (Double.isInfinite(es) || Double.isNaN(es)) {
      Diagnostic.developerLog("s trouble: " + taxonId + " s: " + v + " r: " + r);
      low = 0.0;
      high = 1.0;
    } else {
      low = wt / es;
      final double highZero = wt * es;
      if (highZero > 1.0) {
        Diagnostic.developerLog("s too big: " + taxonId + " s: " + v + " highZero: " + highZero);
        high = 1.0;
      } else {
        high = highZero;
      }
    }
    return new double[] {low, high};
  }

  private double getCoverageDepth(final int taxonId) {
    final SpeciesInfo si = mSpeciesInfo.get(taxonId);
    if (si == null) {
      return 0; // no reads hit this species
    }
    return si.getCoverageDepth();
  }

  private double getCoverageBreadth(final int taxonId) {
    final SpeciesInfo si = mSpeciesInfo.get(taxonId);
    if (si == null) {
      return 0; // no reads hit this species
    }
    return si.getCoverageBreadth();
  }

  private double getMappedReads(final int taxonId) {
    final SpeciesInfo si = mSpeciesInfo.get(taxonId);
    if (si == null) {
      return 0; // no reads hit this species
    }
    return si.getMappedReads();
  }

  void result(SequencesReader templateReader, LineWriter out, BlockResult result, BlockInfo blockInfo) throws IOException {
    final int[] allRef = blockInfo.getSpeciesMap().taxonIds();
    final int[] parentList = makeParentLookup();
    final Vector r = result.getR();
    double tot = 0.0;
    double totDna = 0.0;
    final double[] taxRs = new double[allRef.length];
    final double[] taxDnas = new double[allRef.length];
    final double[] taxMappedReads = new double[allRef.length];
    final long[] taxTotalGenomeLength = new long[allRef.length];
    final long[] taxHitCount = new long[allRef.length];
    final Vector likelihoods = result.getLikelihoods();
    for (int i = 0; i < allRef.length; i++) {
      final double rv = r.get(i);
      final long gl = blockInfo.getGenomeLength(i);
      if (rv > 0.0 || mParams.printAll()) {
        final double dna = rv * gl;
        tot += rv;
        totDna += dna;
        final int taxonId = blockInfo.getTaxonId(i);
        final boolean meetsOutputThreshold = rv > 0 && confidenceDeviations(likelihoods, i) >= parameters().minConfidence();
        final double mappedReads = getMappedReads(taxonId);
        for (int j = i; j >= 0; j = parentList[j]) {
          taxRs[j] += rv;
          taxDnas[j] += dna;
          taxMappedReads[j] += mappedReads;
          taxTotalGenomeLength[j] += gl;
          if (meetsOutputThreshold) {
            taxHitCount[j]++;
          }
        }
      }
    }
    final TreeSet<ComparablePair<ReverseDouble, Integer>> taxset = new TreeSet<>();
    for (int i = 0; i < allRef.length; i++) {
      final double rv = taxRs[i];
      if (rv > 0.0 || mParams.printAll()) {
        taxset.add(new ComparablePair<>(new ReverseDouble(rv), i));
      }
    }
    final Vector varianceLog = result.getVarianceLog();
    double shannon = 0;
    double simpson = 0;

    out.writeln(COMMENT_CHAR + VERSION);
    if (CommandLine.getCommandLine() != null) {
      out.writeln(COMMENT_CHAR + "CL\t" + CommandLine.getCommandLine());
    }
    if (templateReader.getSdfId().available()) {
      out.writeln(COMMENT_CHAR + TEMPLATE_SDF_ATTRIBUTE + TAB + templateReader.getSdfId().toString());
    }

    out.writeln(COMMENT_CHAR + SPECIES_HEADER);
    int speciesCount = 0;
    final HashSet<Integer> tax = mSequenceMap == null ? null : new HashSet<>(mSequenceMap.values());
    final Map<Integer, KronaSpeciesNode> taxonIdToKronaNodeMap = new HashMap<>();
    taxonIdToKronaNodeMap.put(Taxonomy.ROOT_ID, new KronaSpeciesNode(new BoundedDouble(1.0, 1.0, 1.0), new BoundedDouble(1.0, 1.0, 1.0), null, null, null, null, null));
    for (final ComparablePair<ReverseDouble, Integer> co : taxset) {
      final int g = co.getB();
      final double confidence = confidenceDeviations(likelihoods, g);
      final boolean passesConfidenceFilter = confidence >= parameters().minConfidence();
      if (!passesConfidenceFilter && !mParams.printAll()) {
        continue;
      }
      final int taxonId = blockInfo.getTaxonId(g);
      final TaxonNode taxonNode = mTaxonomy.get(taxonId);

      final double w = co.getA().doubleValue();

      final double dna = taxDnas[g];
      final double wtDna = dna / totDna * mSampleFactor;

      final double wt = w / tot;
      // Only the species with genomes will appear in the map
      final boolean hasAssociatedGenome = tax == null || tax.contains(taxonId);
      // Only include the non-printAll species with genomes in the diversity metrics
      // TODO: Actually this might slightly overcount the richness according to one intepretation
      // Any internal node with sequence will be counted if it passes the confidence filter,
      // which can happen even if the internal node is not directly seen, but gets it confidence
      // from children nodes which have been seen.
      if (passesConfidenceFilter && hasAssociatedGenome) {
        shannon -= wt * Math.log(wt);
        simpson += wt * wt;
        speciesCount++;
      }

      final BoundedDouble abundance = outputStdDevBounds(varianceLog, w, g, taxonId, wt);
      final BoundedDouble dnafraction = outputStdDevBounds(varianceLog, w, g, taxonId, wtDna);

      final double coverageDepth = hasAssociatedGenome ? getCoverageDepth(taxonId) / taxTotalGenomeLength[g] : 0;
      final double coverageBreadth = hasAssociatedGenome ? getCoverageBreadth(taxonId) / taxTotalGenomeLength[g] : 0;
      final long genomeLength = hasAssociatedGenome ? taxTotalGenomeLength[g] : 0;

      out.writeln(createSpeciesLine(taxonNode, abundance, dnafraction, confidence, coverageDepth, coverageBreadth, genomeLength, taxMappedReads[g], hasAssociatedGenome, taxHitCount[g]));
      taxonIdToKronaNodeMap.put(taxonId, new KronaSpeciesNode(abundance, dnafraction, confidence, taxMappedReads[g], coverageDepth, coverageBreadth, genomeLength));
    }

    mStatistics.mRichness = speciesCount;
    mStatistics.mShannon = shannon;
    mStatistics.mPielou = shannon / Math.log(speciesCount);
    mStatistics.mInvSimpson = 1.0 / simpson;

    //write Krona HTML stuff if an output dir has been specified (some tests don't)
    //TODO move into statistics report generation
    if (mParams.outputParams() != null && mParams.outputParams().directory() != null) {
      final KronaSpeciesReportWriter ksr = new KronaSpeciesReportWriter(new HtmlReportHelper(mParams.outputParams().directory(), "index"));
      ksr.setParams(mParams);
      ksr.writeReport(mTaxonomy.get(Taxonomy.ROOT_ID), taxonIdToKronaNodeMap);
    }
  }

  private String createSpeciesLine(TaxonNode taxonNode, BoundedDouble abundance, BoundedDouble dnafraction, double confidence, double coverageDepth, double coverageBreadth, long genomeLength, double mappedReads, boolean hasAssociatedGenome, long hitCount) {
    final StringBuilder line = new StringBuilder();

    if (abundance != null) {
      line.append(String.format(Locale.ROOT, NUM_FORMAT, abundance.getValue()));
      line.append(TAB).append(String.format(Locale.ROOT, NUM_FORMAT, abundance.getLow()));
      line.append(TAB).append(String.format(Locale.ROOT, NUM_FORMAT, abundance.getHigh()));
    }
    if (dnafraction != null) {
      line.append(TAB).append(String.format(Locale.ROOT, NUM_FORMAT, dnafraction.getValue()));
      line.append(TAB).append(String.format(Locale.ROOT, NUM_FORMAT, dnafraction.getLow()));
      line.append(TAB).append(String.format(Locale.ROOT, NUM_FORMAT, dnafraction.getHigh()));
    }

    line.append(TAB).append(String.format(Locale.ROOT, NUM_FORMAT_CONFIDENCE, confidence)); // Likelihood AKA confidence

    line.append(TAB).append(String.format(Locale.ROOT, NUM_FORMAT, coverageDepth)); //coverage-depth
    line.append(TAB).append(String.format(Locale.ROOT, NUM_FORMAT, coverageBreadth)); //coverage-breadth
    line.append(TAB).append(genomeLength); //length of reference genome

    line.append(TAB).append(Utils.realFormat(mappedReads, 2)); //mapped-reads
    line.append(TAB).append(hasAssociatedGenome ? "Y" : "N");
    line.append(TAB).append(hitCount);  // number of taxonomy items below and including this node that were hit

    line.append(TAB).append(taxonNode.getId()) //taxonomy-id
            .append(TAB).append(taxonNode.getParentId()) //parent id
            .append(TAB).append(taxonNode.getRank()) //rank
            .append(TAB).append(taxonNode.getName()); // Identifiers
    return line.toString();
  }

  private double confidenceDeviations(final Vector likelihoods, final int g) {
    return Math.sqrt(2.0 * Math.max(0, likelihoods.get(g)));
  }

  private BoundedDouble outputStdDevBounds(Vector varianceLog, double w, int g, int taxonId, double wt) {
    if (varianceLog != null) {
      final double[] bounds = stdDevBounds(taxonId, wt, varianceLog.get(g), w);
      return new BoundedDouble(wt, bounds[0], bounds[1]);
    }
    return null;
  }

  private int[] makeParentLookup() {
    final int[] parentList = new int[mSpeciesMap.size()]; // Each element gives the internal id of the parent
    for (int j = 0; j < parentList.length; j++) {
      final TaxonNode taxon = mTaxonomy.get(mSpeciesMap.taxonId(j));
      final int parentId = taxon.getParentId();
      parentList[j] = parentId != Taxonomy.ROOT_ID ? mSpeciesMap.get(parentId) : -1;
    }
    return parentList;
  }
}
