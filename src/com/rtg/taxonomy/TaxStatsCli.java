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

package com.rtg.taxonomy;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.reader.NamesInterface;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.Counter;
import com.rtg.util.MultiSet;
import com.rtg.util.TextTable;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

/**
 * Taxonomy Verification for SDF
 */
public final class TaxStatsCli extends AbstractCli {

  private static final String SHOW_DETAILS_FLAG = "show-details";

  private final MultiSet<String> mRankCounts = new MultiSet<>(new TreeMap<String, Counter>());
  private final MultiSet<String> mRankLeafCounts = new MultiSet<>(new TreeMap<String, Counter>());
  private Taxonomy mTaxonomy = null;
  private int mLookupNotInTaxonomyCount;
  private int mLookupNoRankCount;
  private int mLookupInternalCount;
  private int mLookupIdsCount;
  private int mLookupSeqsCount;
  private boolean mShowDetails;

  private static class CliValidator implements Validator {
    @Override
    public boolean isValid(CFlags flags) {
      return CommonFlags.validateSDF((File) flags.getAnonymousValue(0));
    }
  }

  @Override
  public String moduleName() {
    return "taxstats";
  }

  @Override
  public String description() {
    return "summarize and verify the taxonomy in an SDF";
  }

  @Override
  protected void initFlags() {
    CommonFlagCategories.setCategories(mFlags);
    mFlags.setDescription("Summarize and perform a verification of a taxonomy within an SDF.");
    mFlags.registerRequired(File.class, CommonFlags.SDF, "SDF to verify taxonomy for").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerOptional(SHOW_DETAILS_FLAG, "list details of sequences attached to internal nodes of the taxonomy").setCategory(CommonFlagCategories.REPORTING);
    mFlags.setValidator(new CliValidator());
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    mShowDetails = mFlags.isSet(SHOW_DETAILS_FLAG);
    return testSdf((File) mFlags.getAnonymousValue(0), new PrintStream(out), err);
  }

  int testSdf(File sdfFile, PrintStream out, PrintStream err) throws IOException {
    try (SequencesReader reader = SequencesReaderFactory.createDefaultSequencesReader(sdfFile)) {
      if (!TaxonomyUtils.hasTaxonomyInfo(reader)) {
        throw new NoTalkbackSlimException("SDF does not contain taxonomy information");
      }
      final Map<String, Integer> sequenceLookupMap = TaxonomyUtils.loadTaxonomyMapping(reader);
      final Set<Integer> lookupIds = new HashSet<>();
      lookupIds.addAll(sequenceLookupMap.values());
      mLookupIdsCount = lookupIds.size();
      mLookupSeqsCount = sequenceLookupMap.keySet().size();
      final Set<String> lookupSequences = sequenceLookupMap.keySet();
      if (!reader.hasNames()) {
        throw new NoTalkbackSlimException("SDF does not have sequence names");
      }
      final NamesInterface names = reader.names();
      final long numSequences = reader.numberSequences();
      if (numSequences > Integer.MAX_VALUE) {
        throw new UnsupportedOperationException();
      }
      final Set<String> sdfSequences = new HashSet<>((int) numSequences);
      for (long i = 0; i < numSequences; ++i) {
        sdfSequences.add(names.name(i));
      }
      mTaxonomy = TaxonomyUtils.loadTaxonomy(reader);
      if (!mTaxonomy.isConsistent()) {
        err.println("SDF taxonomy is inconsistent");
        err.println(mTaxonomy.getInconsistencyReason());
        return 1;
      }
      // Scan for counts and report problems
      int internalSequences = 0;
      int internalNoRank = 0;
      int leafNoRank = 0;
      int leafNodesWithNoSequences = 0;
      // Scan taxonomy collecing counts
      final Map<Integer, Integer> notLeafNodeIds = new HashMap<>();
      for (final TaxonNode current : mTaxonomy.getRoot().depthFirstTraversal()) {
        final String rank = current.getRank();
        final boolean isNoRank = "no rank".equals(rank);
        if (current.isLeaf() && !lookupIds.contains(current.getId())) {
          ++leafNodesWithNoSequences;
        }
        if (!current.isLeaf() && lookupIds.contains(current.getId())) {
          ++internalSequences;
          notLeafNodeIds.put(current.getId(), current.numChildren());
        }
        if (isNoRank) {
          if (current.isLeaf()) {
            ++leafNoRank;
          } else {
            ++internalNoRank;
          }
        }
        mRankCounts.add(rank);
        if (current.isLeaf()) {
          mRankLeafCounts.add(rank);
        }
      }
      // Scan taxon ids mentioned in lookup checking against taxonomy
      for (final Integer id : lookupIds) {
        final TaxonNode current = mTaxonomy.get(id);
        if (current == null) {
          ++mLookupNotInTaxonomyCount;
        } else {
          if (!current.isLeaf()) {
            ++mLookupInternalCount;
          }
          if ("no rank".equals(current.getRank())) {
            ++mLookupNoRankCount;
          }
        }
      }
      // Scan sequence names mentioned in lookup, checking against SDF sequence names
      int numberSequencesMissing = 0;
      final List<String> seqs = new ArrayList<>();
      for (final String lookup : lookupSequences) {
        if (!sdfSequences.contains(lookup)) {
          seqs.add(lookup);
          ++numberSequencesMissing;
        }
      }
      printList(seqs, " sequences in the taxonomy lookup file are not in the SDF", "Sequences in the taxonomy lookup file that are not in SDF:", err);
      // Scan sequence names mentioned in SDF, checking against names in lookup
      seqs.clear();
      for (final String seq : sdfSequences) {
        if (!lookupSequences.contains(seq)) {
          seqs.add(seq);
          ++numberSequencesMissing;
        }
      }
      printList(seqs, " sequences in the SDF are not in the taxonomy lookup file", "Sequences in the SDF that are not in taxonomy lookup file:", err);
      // Inconsistent / incomplete files
      if (numberSequencesMissing > 0) {
        err.println("Error: " + numberSequencesMissing + " sequence names are not in both the SDF and the taxonomy lookup");
        return 1;
      }
      if (mLookupNotInTaxonomyCount > 0) {
        err.println("Error: " + mLookupNotInTaxonomyCount + " taxon ids in the taxonomy lookup are not in the taxonomy");
        return 1;
      }
      // Warnings reporting
      if (internalSequences > 0) {
        err.println("Warning: " + internalSequences + " internal nodes have sequences attached");
        if (mShowDetails) {
          err.println("#nodeId\tchild-count\trank\tname");
          for (final Map.Entry<Integer, Integer> id : notLeafNodeIds.entrySet()) {
            final TaxonNode node = mTaxonomy.get(id.getKey());
            err.println(id.getKey() + "\t" + id.getValue() + "\t" + node.getRank() + "\t" + node.getName());
          }
          err.println();
        }
      }
      if (internalNoRank + leafNoRank > 0) {
        err.println("Warning: " + (internalNoRank + leafNoRank) + " nodes have no rank");
        if (internalNoRank > 0) {
          err.println(internalNoRank + " nodes with no rank are internal nodes");
        }
        if (leafNoRank > 0) {
          err.println(leafNoRank + " nodes with no rank are leaf nodes");
        }
        if (mLookupNoRankCount > 0) {
          err.println(mLookupNoRankCount + " nodes with no rank have sequences attached");
        }
      }
      if (leafNodesWithNoSequences > 0) {
        err.println("Warning: " + leafNodesWithNoSequences + " leaf nodes have no sequences in the SDF");
      }

      // Report general stats
      printStats(out);
      return 0;
    }
  }

  private void printList(List<String> strs, String summaryLabel, String detailLabel, PrintStream err) {
    if (strs.size() > 0) {
      err.println(strs.size() + summaryLabel);
      if (mShowDetails) {
        err.println(detailLabel);
        for (final String str : strs) {
          err.println(str);
        }
        err.println();
      }
    }
  }

  private void printStats(PrintStream out) {
    out.println("TREE STATS");
    final TextTable ttable = new TextTable(1, 0, TextTable.Align.RIGHT);
    ttable.setAlignment(TextTable.Align.LEFT);
    ttable.addRow("internal nodes:", String.valueOf(mTaxonomy.size() - mRankLeafCounts.totalCount()));
    ttable.addRow("leaf nodes:", String.valueOf(mRankLeafCounts.totalCount()));
    ttable.addRow("total nodes:", String.valueOf(mTaxonomy.size()));
    out.println(ttable.toString());
    out.println("RANK COUNTS");
    final TextTable rtable = new TextTable(1, 0, TextTable.Align.RIGHT);
    rtable.setAlignment(TextTable.Align.LEFT);
    rtable.addRow("rank", "internal", "leaf", "total");
    int total = 0;
    int totalLeaf = 0;
    for (final String rank : mRankCounts.keySet()) {
      final int count = mRankCounts.get(rank);
      total += count;
      final int lCount = mRankLeafCounts.get(rank);
      totalLeaf += lCount;
      rtable.addRow("" + rank, String.valueOf(count - lCount), String.valueOf(lCount), "" + count);
    }
    rtable.addRow("TOTAL", String.valueOf(total - totalLeaf), String.valueOf(totalLeaf), String.valueOf(total));
    out.println(rtable.toString());
    out.println("SEQUENCE LOOKUP STATS");
    final TextTable stable = new TextTable(1, 0, TextTable.Align.RIGHT);
    stable.setAlignment(TextTable.Align.LEFT);
    stable.addRow("total sequences:", String.valueOf(mLookupSeqsCount));
    stable.addRow("unique taxon ids:", String.valueOf(mLookupIdsCount));
    stable.addRow("taxon ids in taxonomy:", String.valueOf(mLookupIdsCount - mLookupNotInTaxonomyCount));
    stable.addRow("taxon ids not in taxonomy:", String.valueOf(mLookupNotInTaxonomyCount));
    stable.addRow("internal nodes:", String.valueOf(mLookupInternalCount));
    stable.addRow("leaf nodes:", String.valueOf(mLookupIdsCount - mLookupNotInTaxonomyCount - mLookupInternalCount));
    stable.addRow("no rank nodes:", String.valueOf(mLookupNoRankCount));
    out.println(stable.toString());
  }


}
